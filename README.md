# Barcoded_invasion

Code developed for the DiMartino barcoded invasion project.

## fastaq_read_counter.py

Counts the occurrences of WITS tags in fastaq files.

Only one tag is counted per read.

usage:
$ python3 fastaq_read_counter.py FASTQ_DIR/ WITS_TAG_FILE.fasta OUT_FILE.csv MAX_NO_ERRORS

arg1: FASTQ_DIR/
        Directory that holds fastq files.

arg2: WITS_TAG_FILE.fasta
        FASTA formatted file of tags to count.

arg3: OUT_FILE.csv
        Filename of csv-formatted output file.

arg4: MAX_NO_ERRORS
        (int) maximum number of errors allowed between tag and read

Dependencies: Python (3.5+), Biopython (1.69+), regex (2018.01.10+)

## barcoded_invasion_simulation.py

Simulates #barcoded infection _in silico_