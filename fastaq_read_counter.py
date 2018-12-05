# -*- coding: utf-8 -*-
"""
Counts the occurrences of WITS tags in fastaq files.

Please note that only one tag is counted per read.

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
        

Created on Tue Jun 26 18:55:08 2018

Copyright (c) 2018 Jens Eriksson
"""
import sys
import os
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re
import regex
import csv
import time


def compile_fuzzy_regex_dict(wits_tags, max_n_errors):
    """
    Creates a dictionary with fuzzy regex:es recognizing WITS tag sequences
    
    wits_tags: (list) of fasta seqObjects
    max_n_errors: (int) maximum allowed indels and errors
    
    Returns: (dict) {[seq.id]:(fuzzy regex with n errors allowed)}
    """
    
    fuzzy_regex_dict = {}
    
    for tag in wits_tags:
        regex_to_compile = '('+str(tag.seq.reverse_complement())+'){e<='+str(max_n_errors)+'}'
        fuzzy_regex_dict[tag.id]= regex.compile(regex_to_compile)
    
    return fuzzy_regex_dict

def compile_regular_regex_dict(wits_tags):
    """
    Creates a dict with regular expressions recognizing exact WITS tag sequences. 
    
    RE:s are compiled to the reverse complement of the tag DNA sequence. 
    The RE is stored under a key that is the fasta name (seq.id) of the entries in WITS_tags.
    
    wits_tags: (list) of SeqRecord-objects (read from .fasta file)

    returns: (dict) {[seq.id]:(regex recognizing seq)}
    """
    
    #to store the reverse complement sequence of the WITS tags as regular expressions
    wits_rc_motifs = {}

    for tag in wits_tags:
        regexp = re.compile(str(tag.seq.reverse_complement()))
        wits_rc_motifs[tag.id]=regexp
        
    return wits_rc_motifs
    
def write_csv(fname, counts_dict):
    """
    Writes .csv file that summarizes the counts in counts_dict
    
    fname: (str) name of csv-file
    counts_dict: dict of dicts returned from count_reads()
    """
    
    filenames = list(counts.keys())
    fieldnames = list(wits_rc_motifs.keys()) +  ["not_grouped"] + ["total_reads"] + ["grouped_reads"] + ["exact_matches"] + ["sample"] + ["time_to_process"]
    fieldnames.insert(0, "file")
    with open(fname, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for file in filenames:
            writer.writerow({**{'file':file}, **counts[file]}) #merge dictionaries
    


def count_reads(file_list, data_dir, exact_regex, fuzzy_regex):

    #counts[file]={}, a dict of dicts to store counts for each file
    counts = {} 
    t0 = time.time()

    #Lets go over all the fastq files
    for file in file_list:
        #sample/chip id is the first 10 characters of filename
        sample = file[0:10]
        
        t1 = time.time()
        print("Working on: {}, total time elapsed: {:.1f}".format(file, t1-t0))
        
        #open the file in a proper way
        with open(os.path.join(data_dir, file)) as in_handle:
            
            total = 0 #total reads/sequences in the file
            calls = 0 #assigned tags 
            exact_match = 0 #no exact matches
            counts[file]={} #dict to store total, calls & other calculated values
            
            #iterate through each entry in the loaded .fastq-file
            for title, seq, qual in FastqGeneralIterator(in_handle):
                #store the string representation of the sequence entry
                s=str(seq)
                total += 1 #increment total reads
                
                #print report every 5000 reads
                if (total%5000 == 0):
                    print("Progress...total reads: %i, matches: %i, fails: %i, exact_match: %i" % (total, calls, (total-calls), exact_match))
                
                for wit in exact_regex.keys():
                
                    try_fuzzy_flag = True
                    #test each regex in the exact_regex dict against s
                    if exact_regex[wit].search(s):
                        #if at least one one instance of the tag is in the sequence
                        counts[file][wit] = counts[file].get(wit,0) +1
                        calls += 1
                        exact_match += 1
                        try_fuzzy_flag = False
                        break #no need to test further when match is found
                    
                    if try_fuzzy_flag:
                        #if no exact matches are fount a fuzzy search is initiated
                        if fuzzy_regex[wit].search(s):
                            counts[file][wit] = counts[file].get(wit,0) +1
                            calls +=1
                            break #no need to test further when match is found
                    
            counts[file]["not_grouped"] = total-calls
            counts[file]["total_reads"] = total
            counts[file]["grouped_reads"] = calls
            counts[file]["exact_matches"] = exact_match
            counts[file]["sample"] = sample
            counts[file]["time_to_process_file"] = time.time()-t1
            
            
            print("{} done!, time elapsed for file: {:.1f}".format(file, time.time()-t1))
            
            print("total reads: %i, matches: %i, fails: %i, exact_match: %i" % (total, calls, (total-calls), exact_match))
    
    return counts


def run():
    indir = sys.argv[1]
    wits_file = sys.argv[2]
    csv_name = sys.argv[3]
    n_errors = sys.argv[4]
    
    #working path for fastq-formatted .bam files previously converted with bam2fasta_hack.py
    raw_data_dir = os.path.abspath(indir)
    raw_data_files = os.listdir(raw_data_dir)

    #Load the WITS tags from a FASTA-file in to a list of SeqRecord-objects
    wits_tags = list(SeqIO.parse(wits_file, "fasta"))
    
    #Dict that stores regular (exact) regex
    wits_rc_motifs = compile_regular_regex_dict(wits_tags)
    
    #Dict that stores fuzzy regex 
    fuzzy_wits = compile_fuzzy_regex_dict(wits_tags, n_errors)
    
    results = count_reads(raw_data_files, raw_data_dir, wits_rc_motifs, fuzzy_wits)
    
    write_csv(csv_name, results)

    
    
if __name__ == "__main__"
    run()

        
