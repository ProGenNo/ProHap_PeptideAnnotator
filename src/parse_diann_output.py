import pandas as pd
import argparse
from common import read_fasta
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser(
        description='Reads the DIA-NN report file, parses it into the imput file format suitable for peptide annotation.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input DIA-NN report file")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="full contactenated fasta file")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")       

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()

print ("Reading", args.fasta_file)
fasta_entries = read_fasta(args.fasta_file)

print ("Reading", args.input_file)
diann_df = pd.read_table(args.input_file)
psm_count = len(diann_df)

unique_peptides = diann_df['Stripped.Sequence'].drop_duplicates().tolist()

def find_peptide(seq):
    result_prots = []
    result_pos = []
    for prot in fasta_entries.values():
        if (seq in prot['sequence']):
            result_prots.append(prot['accession'])
            result_pos.append(str(prot['sequence'].index(seq)))

    return [seq, ';'.join(result_prots), ';'.join(result_pos)]

with Pool(args.threads) as p:
    
    pep_map_data = list(tqdm(p.imap_unordered(find_peptide, unique_peptides), total=len(unique_peptides)))
    pep_map = pd.DataFrame(data=pep_map_data, columns=['seq', 'proteins', 'positions']).set_index('seq')

    p.close()
    p.join()

    result_data = []

    for index,row in diann_df.iterrows():
        pept_mapped = pep_map.loc[row['Stripped.Sequence']]
        result_data.append(['row_' + str(index), row['Stripped.Sequence'], pept_mapped['proteins'], pept_mapped['positions']])

    result_df = pd.DataFrame(data=result_data, columns=['ID', 'Sequence', 'Proteins', 'Positions'])
    result_df.to_csv(args.output_file, header=True, index=False, sep='\t')
