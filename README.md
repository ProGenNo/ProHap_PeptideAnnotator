# ProHap Peptide Annotator

Pipeline to annotate peptide identifications when using protein databases made by [ProHap / ProVar](https://github.com/ProGenNo/ProHap).

## Requirements and usage

The pipeline requires Snakemake and Conda installed. You can install these following [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), _Installation via Conda/Mamba_. 

### Steps to annotate a list of PSMs:

1. Format your PSMs or peptides in a tab-separated file as described below
2. Clone the PeptideAnnotator repository: `git clone https://github.com/ProGenNo/ProHap_PeptideAnnotator.git; cd ProHap_PeptideAnnotator`
3. Create a new configuration file based on the instructions in _config_example.yaml_
4. Provide the path to the configuration file on line 2 of _Snakefile_
5. Activate the Conda environment to run Snakemake: `conda activate snakemake`
6. Test Snakemake with a dry-run: `snakemake -c1 -n -q`
7. Run the Snakemake pipeline to perform the annotation, specifying the number of available CPU cores in the `--cores` parameter. E.g., when using 30 cores, run `snakemake --cores 30 -p --use-conda`

### Example

Sample input files and configuration are provided in this repository. Follow these steps to run the example workflow:

```
# Clone the repository
git clone https://github.com/ProGenNo/ProHap_PeptideAnnotator.git;
cd ProHap_PeptideAnnotator;

# Run Snakemake with the default configuration
snakemake --cores 10 -p --use-conda
```

The sample output can be found in the `sample_output.tsv` file.

## Input format

### Required files:

1. List of PSMs or peptides in a tab-separated file having the following four columns (additional columns do not matter):
    - `ID`: Unique identifier for the PSM / peptide
    - `Sequence`: Amino acid sequence of the peptide
    - `Proteins`: List of protein accessions matching the concatenated FASTA file (e.g., `prot_123ab`), separated by semicolon
    - `Positions`: Positions of the first amino acid within the proteins above \(indexed from 0\), separated by semicolon
2. Haplotype table provided by ProHap (if used)
    - If using one of the publicly available ProHap databases, the haplotype table file is provided as the _F2_ file.
3. Variant table provided by ProVar (if used)
    - If using one of the publicly available ProHap databases, ProVar wasn't used and the variant table is not required.
4. Concatenated FASTA file
    - The final protein database produced by ProHap / ProVar. This can be either in the full or simplified format.
5. If using the simplified FASTA format, provide the fasta headers in the separate file.

### Formatting input files:

The input list of peptides or PSMs should be formatted as below:
```
ID	Sequence	Proteins	Positions
pep_fc539	GYEDGGIHLECRSTGWYPQPQIQWSDAK	prot_8489;prot_2c867;prot_4293a;prot_4e288	155;132;113;113
pep_1e5ccd	NYWGSVRR	prot_1003	632
```
The script in `src/format_input.py` can convert any tab- or comma-separated file into the expected format. The script expects the _pandas_, _argparse_, and _tqdm_ libraries (included in the Conda environment in `condaenv.yaml`), and accepts the following arguments:

* `-i`: Input file name (relative or full path)
* `-sep`: Input file separator (default: tab)
* `-f`: FASTA file created by ProHap (full or simplified format)
* `-o`: Output file name (relative or full path)
* `-t`: Number of threads (default: 5)
* `-sc`: Sequence column name (default: "sequence")
* `-id`: ID column name (default: ignore, use line number as unique identifier)
* `-pc`: Protein accessions column name, expecting values separated by colon or semicolon (defualt: ignore, find peptides in the FASTA entries - can be slow)
* `-pos`: Position of peptides within proteins, expecting values separated by colon or semicolon (defualt: ignore, find peptides in the FASTA entries - can be slow)

A provided example shows a simulated output file of [Percolator](http://percolator.ms/) in `sample_data/sample_percolator_output.tsv`. To format this example file, run the following command: 
```
python src/format_input.py -i sample_data/sample_percolator_output.tsv -f sample_data/sample_proteins.fa -sc peptide -id PSMId -pc proteinIds -t 2 -o <path/to/output_file>
```

## Annotation output

The peptide annotation pipeline produces a tab-separated file containing the following columns:

- `ID`: Identifier matching the input file
- `sequence`: Amino acid sequence of the peptide
- `pep_type1`: Peptide classification based on the influence of genetic variation. The possible values are:
    - _contaminant_: peptide matches a contaminant sequence and should be removed from further analysis
    - _canonical_: peptide matches a canonical sequence in Ensembl, without the influence of any genetic variation
    - _single-variant_: peptide contains the product of the alternative allele for one genetic variant
    - _multi-variant_: peptide contains the product of the alternative allele for a set of two or more genetic variants
    - _frameshift_: peptide contains the consequence of the alternative allele of a frameshift variant, the frameshift may also be located upstream of the peptide
    - _variant-no-ref_: Peptide contains the product of the alternative allele for at least one genetic variant. However, the corresponding canonical sequence has not been found in the Ensembl canonical proteome. The peptide should be checked manually.
    - _canonical-no-ref_: Peptide does not cover a product of an alternative allele. However, the sequence has not been found in the Ensembl canonical proteome. The peptide should be checked manually.
- `pep_type2`: Peptide classification based on its ability to distinguish between protein sequences. The possible values are:
    - _contaminant_: peptide matches a contaminant sequence and should be removed from further analysis
    - _proteoform-specific_: maps uniquely to a single form of a protein (i.e., single splice variant and haplotype)
    - _protein-specific_: maps to multiple sequences, which are all products of the same gene
    - _multi-gene_: maps to the products of different genes
- `covered_changes_peptide`: location of matching amino acid changes within the peptide, delimited by a semicolon for multi-variant peptides
- `covered_changes_protein`: covered amino acid changes located within the different matching protein sequences, delimited by "|". E.g., the value `ENST1:123:P>123:R;129:R>129:L|ENST2:223:P>223:R;229:R>229:L` means that there are two amino acid substitutions identified in the _ENST1_ transcript, or two substitutions in the _ENST2_ transcript.
- `covered_alleles_dna`: alleles of genetic variants encoding this peptide. Reference alleles are denoted as _chromosome:position:REF_, alternative alleles denoted as _chromosome:position:REF>ALT_. In case of a multi-gene peptide, the alleles in the respective matching genes will be delimited by "|". E.g., the value `1:123456:C>A;1:123798:G|12:987321:T>G` would mean that in the first matching gene, one alternative and one reference would encode the given peptide, and in the second matching gene, one alternative allele would encode this peptide. The IDs of the respective genes can be found in the _matching_genes_ column.
- `matching_proteins`: identifiers of the proteins matching this peptide. Canonical proteins are identified by their Ensembl Transcript ID (ENSTxxx), contaminants by their provided name (UniProt name in case of cRAP contaminants), and sequences generated by ProHap / ProVar by the identifiers of the respective haplotype or variant sequences.
- `positions_in_proteins`: positions of the peptide within these protein sequences
- `matching_transcripts`: identifiers of the matching transcripts in Ensembl (ENSTxxx)
- `matching_genes`: identifiers of the matching genes in Ensembl (ENSGxxx)
- `preceding_indel_shift`: only relevant for haplotype sequences: the cumulative length of any possible upstream in-frame insertions or deletions. Useful to align peptides with the reference protein
- `reading_frames`: only relevant if three-frame translation in ProHap / ProVar is applied (disabled by default): reading frames used for the translations of the respective protein sequences
- `expected_maximum_frequency`: combined frequency of the haplotypes encoding this peptide sequence. In the case of multi-gene peptides, consider the following example: If the peptide is encoded either by haplotypes A or B of protein 1, or haplotype C of protein 2, the expected maximum frequency will be _max(A.frequency + B.frequency, C.frequency)_.
