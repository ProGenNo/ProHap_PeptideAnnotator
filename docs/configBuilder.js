function getConfigHTML(ensembl_release, max_cores, haplo_table, var_table, final_FASTA, fasta_header, input_file, ID_col, seq_col, prot_col, pos_col, sep, output_file) {
    return "<p>" +
    "# ---------------- General config ----------------<br>" +
    "ensembl_release: " + ensembl_release + "<br>" +
    "max_cores: " + max_cores + "<br>" +
    "<br># ---------------- File paths ----------------<br>" +
    "haplo_db_table: \"" + haplo_table + "\"<br>" +
    "var_db_table: \"" + var_table + "\"<br>" +
    "full_fasta: \"" + final_FASTA + "\"<br>" + 
    "fasta_header: \"" + fasta_header + "\"<br>" +
    "peptides_file: \"" + input_file + "\"<br>" + 
    "output_file: \"" + output_file + "\"<br>" + 
    "<br># ---------------- Input formatting ----------------<br>" +
    "ID_col: \"" + ID_col + "\"<br>" +
    "seq_col: \"" + seq_col + "\"<br>" + 
    "prot_col: \"" + prot_col + "\"<br>" + 
    "pos_col: \"" + pos_col + "\"<br>" +
    "sep: \"" + sep + "\"<br>" +
    "<\p>"
}