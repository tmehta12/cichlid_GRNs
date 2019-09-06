

import os
import sys
import csv
import time
import subprocess

from tqdm import tqdm


# user_input = raw_input("Did you edit the hard coded path? (y/n)")
# if user_input == 'y':
#     continue
# else:
#     sys.exit()


subprocess.run(['mkdir matQualOutput_wg'], shell = True, check = True)
subprocess.run(['mkdir matQualOutput_wg/human'], shell = True, check = True)

file_list = []
for lol in os.listdir(sys.argv[1]):
    # print(lol)
    if '.nfs' not in lol:
        file_list.append(lol)

# print(sorted(file_list, key = lambda x: int(x.split('_')[1])))

for lol in sorted(file_list, key = lambda x: int(x.split('_')[1])):
    print(sys.argv[1] + lol) #motif


    for lel in sorted([ name for name in os.listdir(sys.argv[1] + lol) if os.path.isdir(os.path.join(sys.argv[1] + lol, name)) ]):
        print('\t' + os.path.join(sys.argv[1], lol, lel)) #species
        for lal in sorted(os.listdir(os.path.join(sys.argv[1], lol, lel)), key = lambda x: int(x.split('_')[0])):
            print('\t\t' + os.path.join(sys.argv[1], lol, lel, lal)) #jas, cs, cw
            for lul in sorted(os.listdir(os.path.join(sys.argv[1], lol, lel, lal))):
                if os.path.isdir(os.path.join(sys.argv[1], lol, lel, lal, lul)):
                    print('\t\t\t' + os.path.join(sys.argv[1], lol, lel, lal, lul)) #default, matQualPval
                    for lil in sorted(os.listdir(os.path.join(sys.argv[1], lol, lel, lal, lul))):
                        # print('\t\t\t\t' + os.path.join(sys.argv[1], lol, lel, lal, lul, lil))
                        if lil == 'fimo.txt':
                            with open(os.path.join(sys.argv[1], lol, lel, lal, lul, lil)) as data_file, open('matQualOutput_wg/human/' + '_'.join([lal.split('_')[1], lul, lel, 'results.out']), 'a') as outie:
                                for line in csv.reader(data_file, delimiter = '\t'):
                                    if line[0][0] != '#':
                                        # print(line)
                                        outie.write('\t'.join(line) + '\n')
                                        # time.sleep(.5)
