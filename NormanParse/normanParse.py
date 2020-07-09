from sys import argv, exit
from subprocess import call
import glob
import os
from time import sleep
from subprocess import check_output

'''
    Norman-Parse.py is used to produce and submit a number of slurm batch files that will run 
    tab-parse.py on every file that matches the input path (ARGV[1]).
'''

if len(argv) < 2:
	print("usage: normanParse.py {path containing blast files}")
	exit(0)

file_names = argv[1:]
print(file_names)
print("Running tab-parse.py on {} result files.".format(len(file_names)))
jobid_list =[]
counter = 0

for file in file_names:
    counter += 1
    scriptname = '' + file + '-parse.sh'
    f = open(scriptname, 'w')
    f.write('#!/bin/sh\n')
    f.write('#SBATCH -J sldg-parse_' + str(counter) + '\n')
    f.write('#SBATCH -N 1\n')
    #f.write('#SBATCH -n 48\n')
    f.write('#SBATCH -p defq-48core\n')
    f.write('module load python3/anaconda/5.2.0\n')
    f.write('echo \"run' + str(counter) + '\"\n') 
    f.write('time python tab-parse.py ' + file + '\n')
    f.close()
    os.chmod(scriptname, 0o0777)
    print(scriptname)
    sleep(.1)
    # Output the job id list - needed to check against the job queue to wait for jobs to end.
    jobid = check_output(['sbatch', scriptname])
    jobid = [int(id) for id in jobid.split() if id.isdigit()]
    jobid_list.append(str(jobid).strip('[]'))        

print('\n***  Submitted Parse Jobs  ***\nJob ID(s):')
print('\n'.join(map(str, jobid_list)))

