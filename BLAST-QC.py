import sys
import os
import argparse
import xml.etree.ElementTree as ET

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
# |* Description       : This script is designed to quality control BLAST XML results (BLAST -outfmt 5).
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
        self.accession = None   # <Hit_accession>
        self.length = 0         # <Hit_len>
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
        self.num = 0            # <Iteration_iter-num>
        self.def_ = None        # <Iteration_query-def>
        self.length = 0         # <Iteration_query-len>
        self.hits = []          # Holds a list of all hits for that query
        self.top_hits = []      # Holds a list of the best hits that fit user input

    def order_hits(self, _init_):
        if _init_['order'] == 'e':
            self.hits.sort(key=lambda hit: hit.evalue)
        elif _init_['order'] == 'b':
            self.hits.sort(reverse=True, key=lambda hit: hit.bitscore)
        elif _init_['order'] == 'i':
            self.hits.sort(reverse=True, key=lambda hit: hit.p_identity)
        elif _init_['order'] == 't':
            self.hits.sort(reverse=True, key=lambda hit: hit.deflevel)

        # if a range is specified only add hits that are within range and then sort by deflevel.
        if _init_['erange'] is not 0:
            accept_val = self.hits[0].evalue + _init_['erange']
            for i in range(len(self.hits)):
                if self.hits[i].evalue <= accept_val:
                    self.top_hits.append(self.hits[i])
            self.top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)
        elif _init_['brange'] is not 0:
            accept_val = self.hits[0].bitscore - _init_['brange']
            for i in range(len(self.hits)):
                if self.hits[i].bitscore >= accept_val:
                    self.top_hits.append(self.hits[i])
            self.top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)
        elif _init_['irange'] is not 0:
            accept_val = self.hits[0].p_identity - _init_['irange']
            for i in range(len(self.hits)):
                if self.hits[i].p_identity >= accept_val:
                    self.top_hits.append(self.hits[i])
            self.top_hits.sort(reverse=True, key=lambda hit: hit.deflevel)
        else: self.top_hits = self.hits

        # Apply the input filter number unless all matching hits are desired
        if _init_['num_hits'] != 0:
            self.top_hits = self.top_hits[0:_init_['num_hits']]


class Output:
    def __init__(self, filename):
        self.hits = open(filename+'.hits.txt', 'w')
        self.nohits = open(filename+'.nohits.txt', 'w')
        self.header = open(filename+'.hits.header', 'w')
        self.hits.write("query_name\tquery_length\taccession_number\tsubject_length\tsubject_description\tE value"
                        "\tbit score\tframe\tquery_start\tquery_end\thit_start\thit_end\t%_conserved\t%_identity\n")
        self.nohits.write("query_name\tquery_length\tfailure_reason\n")

    # Print resulting data for all top hits requested into the file in tabular format.
    # Formatting and which values are output can be modified here.
    def write_hits(self, query, _init_):
        for i in range(0, len(query.top_hits)):
            self.hits.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}%\t{}%\n'
                            .format(query.def_, query.length,
                            query.top_hits[i].accession, query.top_hits[i].length, query.top_hits[i].def_,
                            query.top_hits[i].evalue, query.top_hits[i].bitscore, query.top_hits[i].query_frame,
                            query.top_hits[i].query_start, query.top_hits[i].query_end, query.top_hits[i].hit_start,
                            query.top_hits[i].hit_end, query.top_hits[i].p_conserved, query.top_hits[i].p_identity))
            self.header.write('{}\n'.format(query.def_))
        for _ in range(len(query.top_hits), len(query.hits)):
            results_out.nohits.write("{}\t{}\tFiltered by number.\n".format(query.def_, query.length))

    def __del__(self):
        if not self.hits.closed:
            self.hits.close()
        if not self.nohits.closed:
            self.nohits.close()
        if not self.header.closed:
            self.header.close()


# Initializes the argument parser and reads command line args to get filenames, filters and thresholds.
# 'argparse' allows for easy creation of a command line interface menu. If modifactions are desired,
# more info can be found at: 'https://docs.python.org/3/library/argparse.html'
def Initialize():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", help="Specifiy the Blast XML results input file.\n(required)",
                                            required=True, type=str)

    parser.add_argument("-o", "--output", help="Specify the output file base name (no extension). "
                                               "Defaults to base name of input file.", type=str)

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
                                                 "This is defined by how many separate lines are present in the Hit "
                                                 "definition '<Hit_def>' of the XML file.\n(Int value)",
                                                type=int, default=-1)

    parser.add_argument("-or", "--order", help="Specify the order of the results. By lowest evalue, highest bitscore, "
                                                "highest percent identity or most detailed definition data."
                                                "\n(default: by evalue- 'e')",
                                                type=str, choices=["e", "b", "i", "t"], default='e')

    parser.add_argument("-er", "--erange", help="Sets a range of acceptable deviation from the lowest evalue hit "
                                                "in which a more detailed definition would be prefered. "
                                                "Must be ordered by evalue.", type=float, default=0)

    parser.add_argument("-br", "--brange", help="Sets a range of acceptable deviation from the highest bitscore hit "
                                                "in which a more detailed definition would be prefered. "
                                                "Must be ordered by bitscore.", type=float, default=0)

    parser.add_argument("-ir", "--irange", help="Sets a range of acceptable deviation from the highest percent identity "
                                                "hit in which a more detailed definition would be prefered. "
                                                "Must be ordered by percent identity.", type=float, default=0)

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    if args.erange != 0 and args.order != 'e' :
        parser.error('erange cannot be used. Must order by evalue if this functionality is desired.'
                     '\nuse \'-h\' or \'--help\' to display help menu.')
    if args.brange != 0 and args.order != 'b':
        parser.error('brange cannot be used. Must order by bitscore if this functionality is desired.'
                     '\nuse \'-h\' or \'--help\' to display help menu.')
    if args.brange != 0 and args.order != 'i':
        parser.error('brange cannot be used. Must order by bitscore if this functionality is desired.'
                     '\nuse \'-h\' or \'--help\' to display help menu.')

    if not os.path.isfile(args.filename):
        parser.error('The file {} does not exist on this path.'.format(args.filename))

    if not args.filename.lower().endswith('.xml'):
        parser.error('input file must be a BLAST results XML file')

    if args.output is None:
        args.output = args.filename[:-4]

    return {'filename': args.filename, 'output': args.output, 'type': args.type, 'order': args.order,
            'num_hits': args.number, 'bitscore': args.bitscore,  'deflevel': args.definition, 'evalue': args.evalue,
            '%identity': args.identity, 'erange': args.erange, 'brange': args.brange, 'irange': args.irange}

# Taking off . . .
_init_ = Initialize()
results_out = Output(_init_['output'])

# ElementTree is used to parse the XML file to locate data and extract the numerical values from the lines.
# All values are extracted for easy editing of the code, not all are used here. If another value is needed simply
# add another 'cur.find('VALUES_XML-TAG').text' in the desired position of code. More info can be found at:
# 'https://docs.python.org/2/library/xml.etree.elementtree.html'
with open(_init_['filename']) as results_in:
    try:
        tree = ET.parse(results_in)
        root = tree.getroot()
    except:
        raise FileNotFoundError('XML file could not be parsed. Check the BLAST results file: {}.'.format(results_in.name))


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
            if _init_['type'] == 'n':
                cur_hit.deflevel = 1 + cur_hit.def_.count(';')
            if _init_['type'] == 'p':
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
                if cur_hit.evalue <= _init_['evalue'] \
                and cur_hit.bitscore >= _init_['bitscore'] \
                and cur_hit.deflevel >= _init_['deflevel'] \
                and cur_hit.p_identity >= _init_['%identity']:
                    cur_query.hits.append(cur_hit)
                else:
                    results_out.nohits.write("{}\t{}\tBelow Threshold(s).\n".format(cur_query.def_, cur_query.length))

        # If this query had hits apply order and write top hits to result files
        if len(cur_query.hits) != 0:
            cur_query.order_hits(_init_)
            results_out.write_hits(cur_query, _init_)
        else:
            results_out.nohits.write("{}\t{}\tNo hits found.\n".format(cur_query.def_, cur_query.length))
# Touchdown . . .
if not results_in.closed:
    results_in.close()
results_out.__del__()
