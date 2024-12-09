import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Join the Annotator output with the remaining columns from the original file')

parser.add_argument("-orig", dest="orig_file", required=True,
                    help="original report file")

parser.add_argument("-sep", dest="input_sep", required=False,
                    help="input file separator (default: tab)", default='\t')

parser.add_argument("-annot", dest="annot_file", required=True,
                    help="annotated report file")

parser.add_argument("-id", dest="id_col", required=False,
                    help="ID column name (default: N/A (use line number))", default=None)

parser.add_argument("-pc", dest="prot_col", required=False,
                    help="proteins column name (default: N/A (infer))", default=None)

parser.add_argument("-pos", dest="pos_col", required=False,
                    help="position column name (default: N/A (infer))", default=None)

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()

df_orig = pd.read_csv(args.orig_file, sep=args.input_sep)
df_annot = pd.read_table(args.annot_file)

cols_to_drop = []
if (args.prot_col is not None):
    cols_to_drop.append(args.prot_col)
if (args.pos_col is not None):
    cols_to_drop.append(args.pos_col)

result_df = df_annot.join(df_orig.drop(cols_to_drop, axis=1).set_index(args.id_col), on='ID', how='inner', lsuffix='_formatted')
result_df.to_csv(args.output_file, sep='\t', index=False)