import glob
from multiprocessing import cpu_count, Manager, Pool
import csv
import sys        

'''
    tab-parse.py parses a single BLAST tabular results file using multiprocessing across all availible cpu cores.
    hits are compared by bitscore, and only the top bitscore hit is returned for every query sequence.
    the name of the input file is supplied as a command line argument (ARGV[1]).
'''

def parseResult(tabFile, q):
    result = {}
    for row in tabFile:
        if row[0] not in result.keys():
            result[row[0]] = row
        else:
            if(float(result[row[0]][6])<float(row[6])):
                #print("found better bitscore")
                result[row[0]] = row
            else:
                pass
                #print("skip")
    q.put(result)
    #print("got result")

def listener(q, name):
    with open(name[:-3]+'csv', 'w') as f:
        while 1:
            result  = q.get()
            if result == 'kill':
                break
            else:
                for key in result.keys():
                    dataRow=str(key)+','+str(result[key][1])+','+str(result[key][2])+','+str(result[key][3])+','+str(result[key][4])+','+str(result[key][5])+','+str(result[key][6])
                    f.write(dataRow+'\n')
                    #print(dataRow+'\n')
                f.flush()

def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i+n]


filenames = glob.glob(sys.argv[1])
for name in filenames:
    with open(name) as blast_in:
        tabList = list(csv.reader(blast_in, delimiter='\t'))
        splitLists = list(chunks(tabList, cpu_count()))
        manager = Manager()
        q = manager.Queue()
        pool = Pool(processes=cpu_count() + 1)

        pool.apply_async(listener, (q, name))

        processes = []
        for l in splitLists:
            process = pool.apply_async(parseResult, (l,q))
            processes.append(process)

        for process in processes:
            process.get()

        q.put('kill')
        pool.close()
        pool.join()

