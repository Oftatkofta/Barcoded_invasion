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

## barcoded_invasion_simulation_multithreaded.py

Simulates #barcoded infection _in silico_

The simulation is about investigating the stochastic loss of strains in the barcoded invasion assay.

The initial goal is to characterize the influence of invasion probability on
the average number of recovered barcodes.

The model makes the following assumptions:

-  n<sub>cells</sub> = 125.000-175.000 (averages to 150.000)
-  n<sub>bacteria</sub> = MOI*n<sub>cells</sub>
-  P<sub>bind</sub>, the probability that a bacterium binds a cell, is set to 1.0 (all bacteria bind a cell, see comment below).
-  Each bacterium binds one random cell from the cell population.
-  Each cell can have multiple bacteria binding to it.
-  Every bacterium has the probability P<sub>invade</sub> to invade the cell it has bound.
-  If one cell is invaded, all noninvading bacteria bound to the cell have P<sub>coinvade</sub> probability to also invade.
-  P<sub>coinvade</sub> is set to 0, no coinvasion occurs in the simulation.
-  Bacteria do not multiply during the simulation.
-  All invading bacteria are counted with a probability of P<sub>recovered</sub>.
-  P<sub>recovered</sub> is set to 1, all tags are counted.
-  The composition of the innoculum varies on each "infection" with an experimentally derived standard deviation.

The different P<sub>invade</sub> values were determined experimentally. Strictly, this number represents the product of P<sub>bind</sub> and P<sub>invade</sub>.

The model is based on two objects, a Cell, and a Bacterium and is an _in silico_ rewrite of an experimental protocol. It is therefore computationally quite wasteful.

