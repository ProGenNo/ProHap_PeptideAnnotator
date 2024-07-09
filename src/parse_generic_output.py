import pandas as pd
import argparse
import re
from common import read_fasta
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(
        description='Reads the generic PSM report file, parses it into the imput file format suitable for peptide annotation.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input report file")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="full contactenated fasta file")

parser.add_argument("-sc", dest="seq_col", required=False,
                    help="sequence column name (default: 'sequence')", default='sequence')

parser.add_argument("-id", dest="id_col", required=False,
                    help="ID column name (default: 'ID')", default='ID')

parser.add_argument("-pc", dest="prot_col", required=False,
                    help="proteins column name (default: 'proteins')", default=None)

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()

print ("Reading", args.fasta_file)
fasta_entries = read_fasta(args.fasta_file)

print ("Reading", args.input_file)
psm_df = pd.read_table(args.input_file)
psm_count = len(psm_df)

# remove PTMs and other characters (e.g., the N-terminal)
psm_df[args.seq_col] = psm_df[args.seq_col].apply(lambda seq: re.sub(r'\[[^]]*\]', '', seq).split('.')[1].replace('n', '').replace('I', 'L'))

unique_peptides = psm_df[args.seq_col].drop_duplicates().tolist()

def find_peptide(idx):
    row = psm_df.iloc[idx]
    seq = row[args.seq_col]
    result_prots = []
    result_pos = []
    if (args.prot_col is not None):
        for prot in re.split(r"[,;]", row[args.prot_col]):
            if (seq in prot['sequence'].replace('I', 'L')):
                result_prots.append(prot['accession'])
                result_pos.append(str(prot['sequence'].replace('I', 'L').index(seq)))
    else:
        for prot in fasta_entries.values():
            if (seq in prot['sequence'].replace('I', 'L')):
                result_prots.append(prot['accession'])
                result_pos.append(str(prot['sequence'].replace('I', 'L').index(seq)))

    return [row[args.id_col], seq, ';'.join(result_prots) if (len(result_prots) > 0) else row[args.prot_col], ';'.join(result_pos) if (len(result_prots) > 0) else -1]

with Pool(args.threads) as p:
    pep_map_data = list(tqdm(p.imap_unordered(find_peptide, unique_peptides), total=len(unique_peptides)))
    pep_map = pd.DataFrame(data=pep_map_data, columns=['seq', 'proteins', 'positions']).set_index('seq')

    p.close()
    p.join()

    result_data = []

    for index,row in psm_df.iterrows():
        seq = row[args.seq_col]
        pept_mapped = pep_map.loc[seq]
        result_data.append(['row_' + str(index), seq, pept_mapped['proteins'], pept_mapped['positions']])

    # result_data = list(tqdm(p.imap_unordered(find_peptide, range(0, len(psm_df))), total=len(psm_df)))
    result_df = pd.DataFrame(data=result_data, columns=['ID', 'Sequence', 'Proteins', 'Positions'])

    p.close()
    p.join()

    result_df.to_csv(args.output_file, header=True, index=False, sep='\t')
