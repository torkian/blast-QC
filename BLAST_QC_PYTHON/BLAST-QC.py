import sys
import os
import argparse
import xml.etree.ElementTree as ET
import csv
from contextlib import contextmanager
from multiprocessing import Manager, Pool, cpu_count
import traceback

# |**********************************************************************
# |* Project           : Norman Lab Python 3 BLAST Quality Control Script
# |*
# |* Program name      : BLAST-QC.py
# |*
# |* Author            : Spencer Hann, Ben Torkian, Sean Norman
# |*
# |* Date created      : 07/18/2019
# |*
# |* Usage             : python BLAST-QC.py {-args}
# |*
# |* Args              : -h {help} -f {infile} -o {outfile} -t {BLAST type}
# |*                     -n {number of hits} -e {evalue threshold} -b {bit-score threshold}
# |*                     -i {% identity threshold} -d {definition threshold} -or {order by}
# |*                     -er {evalue range} -br {bit-score range} -ir {% identity range}
# |*
# |* description       : This script is designed to quality control BLAST XML results (BLAST -outfmt 5).
# |*                     Results will be filtered based on user input in the form of command-line args,
# |*                     and the best N matching hits with all relevant info are output in a tabular format-
# |*                     for import into a spreadsheet program for analysis.
# |*
# |* License:            This program is free software: you can redistribute it and/or modify
# |*                     it under the terms of the MIT License as published by
# |*                     the Open Source Initiative. 'https://opensource.org/licenses/MIT'
# |*
# |* Revision History  :
# |*
# |* |Date:|        |Author:|      |Ref:|    |Revision:|
# |*
# |*
# |**********************************************************************


class Hit:
    def __init__(self):         # Tag in BLAST XML output
        self.num = 0            # <Hit_num>
        self.id = None          # <Hit_id>
        self.def_ = None        # <Hit_def>
        self.deflevel = 0       # Quantifies the level of information in <Hit_def>
        self.accession = None   # <Hit_accession>
        self.length = 0         # <Hit_len>
        self.mismatch = 0       # number of mismatches (tabular)
        self.gapopen = 0        # number of gap openings (tabular)
        self.bitscore = 0.0     # <Hsp_bitscore>
        self.score = 0          # <Hsp_score>
        self.evalue = 0.0       # <Hsp_evalue>
        self.query_start =0     # <Hsp_query-from>
        self.query_end = 0      # <Hsp_query-to>
        self.hit_start = 0      # <Hsp_hit-from>
        self.hit_end = 0        # <Hsp_hit-to>
        self.query_frame = 0    # <Hsp_query-frame>
        self.positive = 0       # <Hsp_positive>
        self.identity = 0       # <Hsp_identity>
        self.align_len = 0      # <Hsp_align-len>
        self.p_identity = 0.0   # 100*(<Hsp_identity>/<Hsp_align-len>)
        self.p_conserved = 0.0  # 100*(<Hsp_positive>/<Hsp_align-len>)


class Query:
    def __init__(self):
        self.id = None          # qseqid (tabular)
        self.num = 0            # <Iteration_iter-num>
        self.def_ = None        # <Iteration_query-def>
        self.length = 0         # <Iteration_query-len>
        self.hits = []          # Holds a list of all hits for that query


class BLASTQC:
    def __init__(self):
        self.clOptions = self.CLI();

        if self.clOptions.fileformat == "XML":
            with open(self.clOptions.output + '.hits.txt','w') as hits:
                hits.write("query_name\tquery_length\taccession_number\tsubject_length\tsubject_description\tE value"
                            "\tbit score\tframe\tquery_start\tquery_end\thit_start\thit_end\t%_conserved\t%_identity\n")
            with open(self.clOptions.output+'.nohits.txt', 'w') as nohits: 
                nohits.write("query_name\n")
            with open(self.clOptions.output+'.hits.header', 'w') as header: 
                header.write("query_name\tsubject_description\n")
            self.parseXML()

        if self.clOptions.fileformat == "tab":
            with open(self.clOptions.output + '.hits.txt','w') as hits:
                hits.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n")
            with open(self.clOptions.output+'.nohits.txt', 'w') as nohits:
                nohits.write("qseqid\n")
            with open(self.clOptions.output+'.hits.header', 'w') as header: 
                header.write("qseqid\tsseqid\n")
            self.parseTab()


    # Initializes the argument parser and reads command line args to get filenames, filters and thresholds.
    # 'argparse' allows for easy creation of a command line interface menu. If modifactions are desired,
    # more info can be found at: 'https://docs.python.org/3/library/argparse.html'
    def CLI(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--filename", help="Specifiy the Blast XML results input file.\n(required)", type=str)

        parser.add_argument("-ff", "--fileformat", help="Specifiy the Blast results file format (Tabular or XML).\n(required)",
                                                choices=["XML", "tab"], required=True, type=str)

        parser.add_argument("-o", "--output", help="Specify the output file base name (no extension). "
                                                "Defaults to base name of input file.", type=str)

        parser.add_argument("-p", "--parallel", help="Set number of threads for parallel processing. Set to 1 if sequential processing is desired. (Defaults to #of CPU cores avalible.)"
                                                "(INT value)", type=int, default=cpu_count())

        parser.add_argument("-t", "--type", help="Specify what type of BLAST you are running\n(Protein or Nucleotide)."
                                                    " (required)", choices=["p", "n"], type=str, required=True)

        parser.add_argument("-n", "--number", help="Specify the number of hits to return per query sequence. "
                                                "Defaults to return all hits that fit input threshold(s).\n(Int value)",
                                                    type=int, default=0)

        parser.add_argument("-e", "--evalue", help="Specify an e-value threshold.\n(Maximum acceptable evalue)"
                                                "(Float value)", type=float, default=float('Inf'))

        parser.add_argument("-b", "--bitscore", help="Specify a bit-score threshold.\n(Minimum acceptable bitscore)"
                                                    "(Float value)", type=float, default=-1)

        parser.add_argument("-i", "--identity", help="Specify a threshold in the percent identity of a hit."
                                                    "(Calculated value. Not to be confused with identity value)"
                                                    "\n(Minimum acceptable percentage) (Float value)",
                                                    type=float, default=-1)

        parser.add_argument("-d", "--definition", help="Specify a threshold in the level of definition provided."
                                                    "This is defined by the number of line separators (titles) are in the Hit "
                                                    "definition '<Hit_def>' of the XML file, or salltitles column of the tabular output"
                                                    " (must enable the salltitles column in the BLAST tabular output using -outfmt \"6 std salltitles\").\n(Int value)",
                                                    type=int, default=-1)

        parser.add_argument("-or", "--order", help="Specify the order of the results. By lowest evalue, highest bitscore, "
                                                    "highest percent identity or most detailed definition data."
                                                    "\n(default: by evalue- 'e') (if ordering by definition with tabular output you must enable the salltitles column in the BLAST tabular output using -outfmt \"6 std salltitles\")",
                                                    type=str, choices=["e", "b", "i", "d"], default='e')

        parser.add_argument("-er", "--erange", help="Sets a range of acceptable deviation from the lowest evalue hit "
                                                    "in which a more detailed definition would be prefered. "
                                                    "Must be ordered by evalue. (must enable the salltitles column in the BLAST tabular output using -outfmt \"6 std salltitles\" if using tabular output from BLAST)", type=float, default=0)

        parser.add_argument("-br", "--brange", help="Sets a range of acceptable deviation from the highest bitscore hit "
                                                    "in which a more detailed definition would be prefered. "
                                                    "Must be ordered by bitscore. (must enable the salltitles column in the BLAST tabular output using -outfmt \"6 std salltitles\" if using tabular output from BLAST)", type=float, default=0)

        parser.add_argument("-ir", "--irange", help="Sets a range of acceptable deviation from the highest percent identity "
                                                    "hit in which a more detailed definition would be prefered. "
                                                    "Must be ordered by percent identity. (must enable the salltitles column in the BLAST tabular output using -outfmt \"6 std salltitles\" if using tabular output from BLAST)", type=float, default=0)

        args = parser.parse_args()

        if args.erange != 0 and args.order != 'e':
            parser.error('erange cannot be used. Must order by evalue if this functionality is desired.'
                        '\nuse \'-h\' or \'--help\' to display help menu.')
        if args.brange != 0 and args.order != 'b':
            parser.error('brange cannot be used. Must order by bitscore if this functionality is desired.'
                        '\nuse \'-h\' or \'--help\' to display help menu.')
        if args.brange != 0 and args.order != 'i':
            parser.error('brange cannot be used. Must order by identity if this functionality is desired.'
                        '\nuse \'-h\' or \'--help\' to display help menu.')

        if args.output == None and args.filename != None:
            args.output = args.filename[:-4]
        elif args.output == None:
            args.output = "BLASTQC.out"

        return args;

    
    def parseXML(self):
        # ElementTree is used to parse the XML file to locate data and extract the numerical values from the lines.
        # All values are extracted for easy editing of the code, not all are used here. If another value is needed simply
        # add another 'cur.find('VALUES_XML-TAG').text' in the desired position of code. More info can be found at:
        # 'https://docs.python.org/2/library/xml.etree.elementtree.html'
        try:
            tree = None
            if self.clOptions.filename != None:
                tree = ET.parse(open(self.clOptions.filename, 'r'))
            else:
                tree = ET.parse(sys.stdin)
            root = tree.getroot()
        except:
            traceback.print_exc()
            print('\n***  XML file could not be parsed. Check the BLAST results file  ***\n')
            sys.exit()

        for query in root.findall('./BlastOutput_iterations/Iteration'):
            cur_query = Query()
            cur_query.num = query.find('Iteration_iter-num').text
            cur_query.def_ = query.find('Iteration_query-def').text
            cur_query.length = query.find('Iteration_query-len').text
            for hit in query.findall('./Iteration_hits/Hit'):
                cur_hit = Hit()
                cur_hit.id = hit.find('Hit_id').text
                cur_hit.def_ = hit.find('Hit_def').text

                # Change formating of <Hit_def> with blast type - different delimiters.
                # deflevel is defined by the count of those delimiters, as with each one there is
                # an increase in the level of detail in the definition of the hit.
                if self.clOptions.type == 'n':
                    cur_hit.deflevel = 1 + cur_hit.def_.count(';')
                if self.clOptions.type == 'p':
                    cur_hit.deflevel = 1 + cur_hit.def_.count('>')
                cur_hit.accession = hit.find('Hit_accession').text
                cur_hit.length = hit.find('Hit_len').text
                count = 1
                
                for hsp in hit.findall('./Hit_hsps/Hsp'):
                    # Count value is needed to solve the case of multiple hsps per hit.
                    # Create a new hit object and treat it as a separate hit in the list although it retains some values.
                    if count > 1:
                        new_hit = Hit()
                        new_hit.id = cur_hit.id
                        new_hit.def_ = cur_hit.def_
                        new_hit.deflevel = cur_hit.deflevel
                        new_hit.accession = cur_hit.accession
                        new_hit.length = cur_hit.length
                        cur_hit = new_hit
                        
                    cur_hit.bitscore = float(hsp.find('Hsp_bit-score').text)
                    cur_hit.score = int(hsp.find('Hsp_score').text)
                    cur_hit.evalue = float(hsp.find('Hsp_evalue').text)
                    cur_hit.query_start =int(hsp.find('Hsp_query-from').text)
                    cur_hit.query_end = int(hsp.find('Hsp_query-to').text)
                    cur_hit.hit_start = int(hsp.find('Hsp_hit-from').text)
                    cur_hit.hit_end = int(hsp.find('Hsp_hit-to').text)
                    cur_hit.query_frame = int(hsp.find('Hsp_query-frame').text)
                    cur_hit.identity = float(hsp.find('Hsp_identity').text)
                    cur_hit.align_len = float(hsp.find('Hsp_align-len').text)
                    cur_hit.positive = float(hsp.find('Hsp_positive').text)
                    count += 1

                    # Calculate the %identity and %conserved by using the align length and identity/positive data
                    cur_hit.p_identity = float("%.1f"%(100 * cur_hit.identity / cur_hit.align_len))
                    cur_hit.p_conserved = float("%.1f"%(100 * cur_hit.positive / cur_hit.align_len))

                    # Apply thresholds and add hit to list if it conforms.
                    # If changes to the thresholds are desired (add more ect.) this is where do do it.
                    if cur_hit.evalue <= self.clOptions.evalue \
                    and cur_hit.bitscore >= self.clOptions.bitscore \
                    and cur_hit.deflevel >= self.clOptions.definition \
                    and cur_hit.p_identity >= self.clOptions.identity:
                        cur_query.hits.append(cur_hit)

            # If this query had hits apply order and write top hits to result files
            if len(cur_query.hits) != 0:
                self.order_hits(cur_query)
                with open(self.clOptions.output+'.hits.txt', 'a') as hits:
                    with open(self.clOptions.output+'.hits.header', 'a') as header:
                        for i in range(0, len(cur_query.hits)):
                            hits.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}%\t{}%\n'
                                            .format(cur_query.def_, cur_query.length,
                                            cur_query.hits[i].accession, cur_query.hits[i].length, cur_query.hits[i].def_,
                                            cur_query.hits[i].evalue, cur_query.hits[i].bitscore, cur_query.hits[i].query_frame,
                                            cur_query.hits[i].query_start, cur_query.hits[i].query_end, cur_query.hits[i].hit_start,
                                            cur_query.hits[i].hit_end, cur_query.hits[i].p_conserved, cur_query.hits[i].p_identity))
                            header.write('{}\t{}\n'.format(cur_query.def_, cur_query.hits[i].def_))
            else:
                with open(self.clOptions.output+'.nohits.txt', 'a') as nohits:
                    nohits.write("{}\tNo hits found.\n".format(cur_query.def_))
            

    def parseTab(self):
        # csv is used to parse the blast tabular output format (outfmt 6)
        # more info may be found at 'https://docs.python.org/3/library/csv.html'.
        # if using parsing options that rely on the definition of a hit, enable the column 'salltitles' in -outfmt 6 of BLAST
        # e.g. '-outfmt "6 std salltitles"'
        if(self.clOptions.filename != None):
            results_in = open(self.clOptions.filename, 'r')
        else:
            results_in = sys.stdin

        cur_query = None;
        cur_hit = None;
        tabfile = csv.reader(results_in, delimiter='\t')
        for row in tabfile:
            if(row[0][0] == '#'):
                continue
            elif cur_query == None:
                cur_query = Query()
            elif cur_query.id != row[0]:
                if len(cur_query.hits) != 0:
                    self.order_hits(cur_query)
                    with open(self.clOptions.output+'.hits.txt', 'a') as hits:
                        with open(self.clOptions.output+'.hits.header', 'a') as header:
                            for i in range(0, len(cur_query.hits)):
                                hits.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
                                                .format(cur_query.id, cur_query.hits[i].id, cur_query.hits[i].identity, 
                                                cur_query.hits[i].align_len, cur_query.hits[i].mismatch, cur_query.hits[i].gapopen, cur_query.hits[i].query_start,
                                                cur_query.hits[i].query_end, cur_query.hits[i].hit_start, cur_query.hits[i].hit_end, cur_query.hits[i].evalue, 
                                                cur_query.hits[i].bitscore, cur_query.hits[i].def_))
                                header.write('{}\t{}\n'.format(cur_query.id, cur_query.hits[i].id))
                else:
                    with open(self.clOptions.output+'.nohits.txt', 'a') as nohits:
                        nohits.write("{}\tNo hits found.\n".format(cur_query.id))
                cur_query = Query()

            cur_hit = Hit()
            cur_query.id = row[0]
            cur_hit.id = row[1]
            cur_hit.identity = row[2]
            cur_hit.align_len = row[3]
            cur_hit.mismatch = row[4]
            cur_hit.gapopen = row[5]
            cur_hit.query_start = row[6]
            cur_hit.query_end = row[7]
            cur_hit.hit_start = row[8]
            cur_hit.hit_end = row[9]
            cur_hit.evalue = row[10]
            cur_hit.bitscore = row[11]

            if len(row) >= 13:
                cur_hit.def_ = row[12]
                if self.clOptions.type == 'n':
                    cur_hit.deflevel = 1 + cur_hit.def_.count(';')
                elif self.clOptions.type == 'p':
                    cur_hit.deflevel = 1 + cur_hit.def_.count('>')

             # Calculate the %identity and %conserved by using the align length and identity/positive data
            cur_hit.p_identity = float("%.1f"%(100 * cur_hit.identity / cur_hit.align_len))

            cur_query.hits.append(cur_hit)


    def order_hits(self, cur_query):
        if(self.clOptions.parallel > 1):
            cur_query.hits = self.parallel_merge_sort(cur_query.hits, self.clOptions.parallel, self.clOptions.order)
        else:
            if self.clOptions.order == 'e':
                cur_query.hits.sort(key=lambda hit: hit.evalue)
            elif self.clOptions.order == 'b':
                cur_query.hits.sort(reverse=True, key=lambda hit: hit.bitscore)
            elif self.clOptions.order == 'i':
                cur_query.hits.sort(reverse=True, key=lambda hit: hit.p_identity)
            elif self.clOptions.order == 'd':
                cur_query.hits.sort(reverse=True, key=lambda hit: hit.deflevel)

        top_hits = [];
        # if a range is specified only add hits that are within range and then sort by deflevel.
        if self.clOptions.erange is not 0:
            accept_val = cur_query.hits[0].evalue + self.clOptions.erange
            for i in range(len(cur_query.hits)):
                if cur_query.hits[i].evalue <= accept_val:
                    top_hits.append(cur_query.hits[i]) 
                else:
                    break  
            if(self.clOptions.parallel > 1):
                top_hits = self.parallel_merge_sort(top_hits, self.clOptions.parallel, 'd')
            else:
                top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)

        elif self.clOptions.brange is not 0:
            accept_val = cur_query.hits[0].bitscore - self.clOptions.brange
            for i in range(len(cur_query.hits)):
                if cur_query.hits[i].bitscore >= accept_val:
                    top_hits.append(self.hits[i])
                else:
                    break 
            if(self.clOptions.parallel > 1):
                top_hits = self.parallel_merge_sort(top_hits, self.clOptions.parallel, 'd')
            else:
                top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)

        elif self.clOptions.irange is not 0:
            accept_val = cur_query.hits[0].p_identity - self.clOptions.irange
            for i in range(len(cur_query.hits)):
                if cur_query.hits[i].p_identity >= accept_val:
                    top_hits.append(cur_query.hits[i])
                else:
                    break
            if(self.clOptions.parallel > 1):
                top_hits = self.parallel_merge_sort(top_hits, self.clOptions.parallel, 'd')
            else:
                top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)

        else: top_hits = cur_query.hits

        # Apply the input filter number unless all matching hits are desired
        if self.clOptions.number != 0:
            cur_query.hits = cur_query.hits[0:self.clOptions.number]


    def comp(self, hit1, hit2, mode):
        if mode == 'e':
            return hit1.evalue < hit2.evalue
        elif mode == 'b':
            return hit1.bitscore > hit2.bitscore
        elif mode == 'i':
            return hit1.p_identity > hit2.p_identity
        elif mode == 'd':
            return hit1.deflevel > hit2.deflevel


    def merge_sort_multiple(self, results, array, mode):
        results.append(self.merge_sort(array, mode))


    def merge_multiple(self, results, array_part_left, array_part_right, mode):
        results.append(self.merge(array_part_left, array_part_right, mode))


    def merge_sort(self, array, mode):
        array_length = len(array)

        if array_length <= 1:
            return array

        middle_index = int(array_length / 2)
        left = array[0:middle_index]
        right = array[middle_index:]
        left = self.merge_sort(left, mode)
        right = self.merge_sort(right, mode)
        return self.merge(left, right, mode)


    def merge(self, left, right, mode):
        sorted_list = []
        # We create shallow copies so that we do not mutate
        # the original objects.
        left = left[:]
        right = right[:]
        # We do not have to actually inspect the length,
        # as empty lists truth value evaluates to False.
        # This is for algorithm demonstration purposes.
        while len(left) > 0 or len(right) > 0:
            if len(left) > 0 and len(right) > 0:
                if self.comp(left[0], right[0], mode):
                    sorted_list.append(left.pop(0))
                else:
                    sorted_list.append(right.pop(0))
            elif len(left) > 0:
                sorted_list.append(left.pop(0))
            elif len(right) > 0:
                sorted_list.append(right.pop(0))
        return sorted_list


    @contextmanager
    def process_pool(self, size):
        """Create a process pool and block until
        all processes have completed.
        Note: see also concurrent.futures.ProcessPoolExecutor"""
        pool = Pool(size)
        yield pool
        pool.close()
        pool.join()


    def parallel_merge_sort(self, array, process_count, mode):
        # Divide the list in chunks
        step = int(len(array) / process_count)

        # Instantiate a multiprocessing.Manager object to
        # store the output of each process.
        # See example here
        # http://docs.python.org/library/multiprocessing.html#sharing-state-between-processes
        manager = Manager()
        results = manager.list()

        with self.process_pool(process_count) as pool:
            for n in range(process_count):
                # We create a new Process object and assign the
                # merge_sort_multiple function to it,
                # using as input a sublist
                if n < process_count - 1:
                    chunk = array[n * step:(n + 1) * step]
                else:
                    # Get the remaining elements in the list
                    chunk = array[n * step:]
                pool.apply_async(self.merge_sort_multiple, (results, chunk, mode))

        # For a core count greater than 2, we can use multiprocessing
        # again to merge sub-lists in parallel.
        while len(results) > 1:
            with self.process_pool(size=process_count) as pool:
                pool.apply_async(
                    self.merge_multiple,
                    (results, results.pop(0), results.pop(0), mode)
                )

        final_sorted_list = results[0]

        return final_sorted_list



# Taking off . . .
BLASTQC()