
# returns an object containing all the proteins + metadata in the fasta file,  
# accessed by the stable protein ID or accession in case of artificial identifier
# all_proteins[proteinID] = {'tag': tag, 'accession': accession, 'description': description, 'sequence': sequence}
# if fasta headers are stored in a separate file, provide the data frame with the headers as the second argument, with the protein accession as index
def read_fasta(filename, headers=None, trim_acc=False):
    proteinDB_file = open(filename, 'r')

    # read protein db
    metadata = proteinDB_file.readline()    # line starting with '>'
    sequence = ""                           # all the other lines are considered protein sequences (considering also a multi-line format)

    all_proteins = {}

    # read the reference sequence database
    while metadata != "":

        line = proteinDB_file.readline()
        while not line.startswith('>') and not line == "":
            sequence += line[:-1]
            line = proteinDB_file.readline()

        tag = ''
        accession = ''
        description = ''

        if "|" in metadata:					# the header is at least partially formated
            metadata_parsed = metadata[1:].split('|')
            if (headers is not None):
                tag = headers.loc[accession]['tag']
            elif 'generic' in metadata_parsed[0]:
                tag = metadata_parsed[0]
            else:
                tag = 'generic_' + metadata_parsed[0]

            if len(metadata_parsed) == 2:					# accession and description are potentially merged -> separate them
                if " " in metadata_parsed[1]:
                    accession = metadata_parsed[1].split(' ')[0]
                    description = metadata_parsed[1].split(' ', 1)[1]
                else:
                    accession = metadata_parsed[1]				# no description -> keep accession as it is
            elif len(metadata_parsed) == 3:					# descripton and accesson already separated -> keep
                accession = metadata_parsed[1]
                description = metadata_parsed[2]

        else:						                    # the header is not formated
            accession = metadata[1:-1].split(" ")[0]
            if " " in metadata:
                description = metadata.split(" ", 1)[1]
            if (headers is not None):
                tag = headers.loc[accession]['tag']

        if trim_acc:
            proteinID = accession.split('.')[0]
        else:
            proteinID = accession

        matching_proteins = []
        seq_positions = []
        reading_frames = []

        if 'position_within_protein:' in description:
            seq_positions = [ int(pos) for pos in description.split('position_within_protein:', 1)[1].split(' ', 1)[0].split(';') ]
        elif (headers is not None):
            seq_positions = [ int(pos) for pos in headers.loc[accession]['position_within_protein'].split(';') ]

        if 'matching_proteins:' in description:
            proteinLists = description.split('matching_proteins:', 1)[1].split(' ', 1)[0].split(';')
            matching_proteins = [ l.split(',') for l in proteinLists ]
        elif (headers is not None):
            matching_proteins = [ l.split(',') for l in headers.loc[accession]['matching_proteins'].split(';') ]

        if 'reading_frame:' in description:
            rfLists = description.split('reading_frame:', 1)[1].split(maxsplit=1)[0].split(';')
            reading_frames = [ l.split(',') for l in rfLists ]
        elif (headers is not None):
            reading_frames = [ l.split(',') for l in headers.loc[accession]['reading_frame'].split(';') ]

        all_proteins[proteinID] = {'tag': tag, 'accession': accession, 'description': description.replace('\n', ''), 'sequence': sequence, 'seq_positions': seq_positions, 'matching_proteins': matching_proteins, 'reading_frames': reading_frames}

        metadata = line
        sequence = ""

    proteinDB_file.close()

    return all_proteins
