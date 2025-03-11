# Specify the configuration file below:
configfile: "config.yaml"

Ensembl_FTP_URL = "ftp.ensembl.org/pub/release-" + str(config['ensembl_release']) + "/"
annotationFilename = "Homo_sapiens.GRCh38." + str(config['ensembl_release']) + ".chr_patch_hapl_scaff"

rule all:
    input:
        config['output_file']

rule download_reference_proteome:
    output:
        temp("data/fasta/Homo_sapiens.GRCh38.pep.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + Ensembl_FTP_URL + "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O {output}.gz && gunzip {output}.gz; "

rule reference_fix_headers:
    input:
        "data/fasta/Homo_sapiens.GRCh38.pep.all.fa"
    output:
        "data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa"
    conda: "condaenv.yaml"
    shell:
        "python3 src/fix_headers.py -i {input} -o {output} -t _ensref -use_ENST 1 "

rule download_gtf:
    output:
        temp("data/gtf/" + annotationFilename + ".gtf")
    shell:
        "mkdir -p data/gtf ; "
        "wget " + Ensembl_FTP_URL + "gtf/homo_sapiens/" + annotationFilename + ".gtf.gz -O {output}.gz && gunzip {output}.gz; "

rule parse_gtf_whole:
    input:
        "data/gtf/" + annotationFilename + ".gtf"
    output:
        "data/gtf/" + annotationFilename + ".db"
    conda: "condaenv.yaml"
    shell:
        "python src/parse_gtf.py -i {input} -o {output}"

rule format_input:
    input:
        pep=config['peptides_file'],
        fasta_file=config['full_fasta'],
    output:
        pep=temp("data/peptides_formatted.tsv"),
        log=".".join(config['output_file'].split('.')[:-1]) + '_PeptideAnnotator_log.txt'
    params:
        id_col=config['ID_col'],
        seq_col=config['seq_col'],
        prot_col=config['prot_col'],
        pos_col=config['pos_col'],
        sep=config['sep'],
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "condaenv.yaml"
    shell:
        "mkdir -p data ; python src/format_input.py -i \"{input.pep}\" -f \"{input.fasta_file}\" -id {params.id_col} -sc {params.seq_col} " +
        ("-sep {params.sep} " if (config['sep'] != "\t") else "") +
        ("-pc {params.prot_col} " if (len(config['prot_col']) > 0) else "") +
        ("-pos {params.pos_col} " if (len(config['pos_col']) > 0) else "") +
        (("-dc \"" + str(config['decoy_col']) + "\" -dv " + str(config['decoy_val']) + " ") if (len(config["decoy_col"]) > 0) else "") +
        "-t {params.max_cores} -removed \"{output.log}\" -o \"{output.pep}\""    

rule annotate_peptides:
    input:
        pep="data/peptides_formatted.tsv",
        var_db=expand('{proxy}', proxy=[config['var_db_table']] if len(config["var_db_table"]) > 0 else []),
        haplo_db=expand('{proxy}', proxy=[config['haplo_db_table']] if len(config["haplo_db_table"]) > 0 else []),
        annot_db="data/gtf/" + annotationFilename + ".db",
        fasta_file=config['full_fasta'],
        fasta_header=expand('{proxy}', proxy=[config['fasta_header']] if len(config["fasta_header"]) > 0 else []),
        ref_fasta='data/fasta/ensembl_reference_proteinDB_' + str(config['ensembl_release']) + '_tagged.fa'
    output:
        temp("data/peptides_annotated.tsv")
    params:
        variant_prefix="var_",
        haplotype_prefix="haplo_",
        max_cores=config['max_cores'],
        log_file=".".join(config['output_file'].split('.')[:-1]) + '_PeptideAnnotator_log.txt',
        header_arg=("-fh \"" + config['fasta_header'] + '\"') if len(config["fasta_header"]) else ""
    threads: config['max_cores']
    conda: "condaenv.yaml"
    shell:
        "python src/peptides_annotate_variation.py -i {input.pep} " +
        ("-var_tsv \"{input.var_db}\" -var_prefix {params.variant_prefix} " if (len(config["var_db_table"]) > 0) else "") +
        ("-hap_tsv \"{input.haplo_db}\" -hap_prefix {params.haplotype_prefix} " if (len(config["haplo_db_table"]) > 0) else "") +
        "-decoy_col \"Decoy\" -decoy_val 1 -log \"{params.log_file}\" -ens_annot \"{input.annot_db}\" -f \"{input.fasta_file}\" {params.header_arg} -ref_fa \"{input.ref_fasta}\" -t {params.max_cores} -o \"{output}\"; "

rule join_output:
    input:
        orig=config['peptides_file'],
        annot="data/peptides_annotated.tsv"
    output:
        config['output_file']
    params:
        id_col=config['ID_col'],
        prot_col=config['prot_col'],
        pos_col=config['pos_col'],
        sep=config['sep']    
    conda: "condaenv.yaml"
    shell:
        "python src/join_reports.py -orig \"{input.orig}\" -annot \"{input.annot}\" -id \"{params.id_col}\" " +
        ("-sep \"{params.sep}\" " if (config['sep'] != "\t") else "") +
        ("-pc \"{params.prot_col}\" " if (len(config['prot_col']) > 0) else "") +
        ("-pos \"{params.pos_col}\" " if (len(config['pos_col']) > 0) else "") +
        "-o \"{output}\""