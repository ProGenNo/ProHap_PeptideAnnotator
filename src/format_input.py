import pandas as pd
import argparse
import re
from datetime import datetime
from common import read_fasta
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(
        description='Reads the generic PSM report file, parses it into the imput file format suitable for peptide annotation.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input report file")

parser.add_argument("-sep", dest="input_sep", required=False,
                    help="input file separator (default: tab)", default='\t')

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="full contactenated fasta file")

parser.add_argument("-sc", dest="seq_col", required=False,
                    help="sequence column name (default: 'sequence')", default='sequence')

parser.add_argument("-id", dest="id_col", required=False,
                    help="ID column name (default: N/A (use line number))", default=None)

parser.add_argument("-pc", dest="prot_col", required=False,
                    help="proteins column name (default: N/A (infer))", default=None)

parser.add_argument("-pos", dest="pos_col", required=False,
                    help="position column name (default: N/A (infer))", default=None)

parser.add_argument("-t", dest="threads", type=int, required=False,
                    help="maximum number of threads (default: 5)", default=5)

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

parser.add_argument("-removed", dest="removed_output_file", required=True,
                    help="output TSV file - list of removed peptides (not matching provided FASTA)")

args = parser.parse_args()

print ("Reading", args.fasta_file)
fasta_entries = read_fasta(args.fasta_file)

print ("Reading", args.input_file)
psm_df = pd.read_csv(args.input_file, sep=args.input_sep)
psm_count = len(psm_df)

# remove PTMs and other characters (e.g., the N-terminal, charge state)
# special condition for Percolator format to remove the residues before and after peptide (e.g., M.n[+42.021]PEPTIDEK2.P -> PEPTIDEK)
psm_df[args.seq_col] = psm_df[args.seq_col].apply(lambda seq: re.sub(r'\[[^]]*\]|\d', '', seq).split('.')[1].replace('n', '').replace('I', 'L') if (seq[1] == '.') else re.sub(r'\[[^]]*\]', '', seq).replace('I', 'L') )

unique_peptides = psm_df.drop_duplicates(subset=[args.seq_col])

def find_peptide(idx):
    row = unique_peptides.iloc[idx]
    seq = row[args.seq_col]
    result_prots = []
    result_pos = []
    if (args.prot_col is not None):
        if (args.pos_col is not None):
            result_pos = re.split(r"[,;]", row[args.pos_col])
            result_prots = row[args.prot_col].split(';')
        else:
            for protID in re.split(r"[,;]", row[args.prot_col]):
                prot = fasta_entries[protID]
                if (seq in prot['sequence'].replace('I', 'L')):
                    result_prots.append(prot['accession'])
                    result_pos.append(str(prot['sequence'].replace('I', 'L').index(seq)))
    else:
        for prot in fasta_entries.values():
            if (seq in prot['sequence'].replace('I', 'L')):
                result_prots.append(prot['accession'])
                result_pos.append(str(prot['sequence'].replace('I', 'L').index(seq)))

    return [seq, ';'.join(result_prots) if (len(result_prots) > 0) else (row[args.prot_col] if (args.prot_col is not None) else '-'), ';'.join(result_pos) if (len(result_prots) > 0) else -1]

with Pool(args.threads) as p:
    pep_map_data = list(tqdm(p.imap_unordered(find_peptide, range(0, len(unique_peptides))), total=len(unique_peptides)))
    pep_map = pd.DataFrame(data=pep_map_data, columns=['seq', 'proteins', 'positions']).set_index('seq')

    p.close()
    p.join()

    result_data = []

    for index,row in psm_df.iterrows():
        seq = row[args.seq_col]
        pept_mapped = pep_map.loc[seq]
        result_data.append([row[args.id_col] if (args.id_col is not None) else 'pep_' + hex(index)[2:], seq, pept_mapped['proteins'], pept_mapped['positions']])

    # result_data = list(tqdm(p.imap_unordered(find_peptide, range(0, len(psm_df))), total=len(psm_df)))
    result_df = pd.DataFrame(data=result_data, columns=['ID', 'Sequence', 'Proteins', 'Positions'])

    removed_outfile = open(args.removed_output_file, 'w')
    removed_outfile.write('------------' + '[' + datetime.now().strftime('%X %x') + '] file: ' + args.input_file + ' ------------\nRemoved peptides:\n')
    for index,row in result_df[result_df['Proteins'] == '-'].iterrows():
        removed_outfile.write(row['ID'] + ': ' + row['Sequence'] + '\n')
    removed_outfile.close()
    
    # Remove peptides that do not match any protein in the provided FASTA
    result_df = result_df[result_df['Proteins'] != '-']

    p.close()
    p.join()

    result_df.to_csv(args.output_file, header=True, index=False, sep='\t')
