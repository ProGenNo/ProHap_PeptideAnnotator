# Specify the configuration file below:
configfile: "config_example.yaml"

Ensembl_FTP_URL = "ftp.ensembl.org/pub/release-" + str(config['ensembl_release']) + "/"
annotationFilename = "Homo_sapiens.GRCh38." + str(config['ensembl_release']) + ".chr_patch_hapl_scaff"

rule all:
    input:
        config['output_file']

rule download_reference_proteome:
    output:
        "data/fasta/Homo_sapiens.GRCh38.pep.all.fa"
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

rule annotate_peptides:
    input:
            pep=config['peptides_file'],
            var_db=expand('{proxy}', proxy=[config['var_db_table']] if len(config["var_db_table"]) > 0 else []),
            haplo_db=expand('{proxy}', proxy=[config['haplo_db_table']] if len(config["haplo_db_table"]) > 0 else []),
            annot_db="data/gtf/" + annotationFilename + ".db",
            fasta_file=config['full_fasta'],
            ref_fasta='data/fasta/ensembl_reference_proteinDB_' + str(config['ensembl_release']) + '_tagged.fa'
    output:
            config['output_file']
    params:
            variant_prefix="var_",
            haplotype_prefix="haplo_",
            max_cores=config['max_cores'],
            log_file="annotation_log.txt"
    threads: config['max_cores']
    conda: "condaenv.yaml"
    shell:
            "python src/peptides_annotate_variation.py -i {input.pep} " +
            ("-var_tsv {input.var_db} -var_prefix {params.variant_prefix} " if (len(config["var_db_table"]) > 0) else "") +
            ("-hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} " if (len(config["haplo_db_table"]) > 0) else "") +
            "-log {params.log_file} -ens_annot {input.annot_db} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "
