#!/usr/bin/python
#=====================================================================================================================================================
'''
Parses results from rsat matrix_quality.pl when run (with 1000 matrix permutations per motif) on a list of transcription factors

Will create a matirx of pvalues for use in fimo scanning

'''
print('\n===================================================================================\n')
import sys

print('matQualResultsParser1.py\n')

if len(sys.argv) < 3:
    print('\nNot enough arguments entered.\nUsage: JASPAR File, hsap JASPAR Endsembl, mmus JASPAR Ensembl, Dir of rsat results dirs\n')
    print('\n===================================================================================\n')
    sys.exit()


print("Importing...\n\n")
import os
import csv
import time
import string
# import Bio
# import numpy
# import itertools
import subprocess
import collections

# print(sys.argv)
# from Bio import AlignIO

# print('!!!\tRequrires biopython & mafft loading\t!!!\n\n')

from tqdm import tqdm
from collections import OrderedDict


####################################################
####################################################
dir_list = []
for dire in [ name for name in os.listdir('RG-cich/CSPSSMs') if os.path.isdir(os.path.join('RG-cich/CSPSSMs', name)) ]:

# for dire in os.listdir('RG-cich/CSPSSMs/'):
    # print(dire)
    if 'mm10' not in dire:
        if 'hE' not in dire:
            if 'mE' not in dire:
                if 'logs' not in dire:
                    if 'sub' not in dire:
                        dir_list.append(dire)
    # if os.path.isdir(dire):
    #     print(dire)
    #
    #     if dire not in cich_motifs_list:
    #         print('\t' + dire)
####################################################
####################################################

def parseJASPAR(JasparTransfacs):
    motif_dict  = OrderedDict()
    with open(JasparTransfacs) as motif_file:
        lines = motif_file.read()
        for motif in lines.split('//'):

            for line in motif.split('\n'):

                if line.startswith('AC'):
                    motif_symbol             = line.replace('AC ', '')

                if line.startswith('ID'):
                    motif_dict[line.replace('ID ', '').upper()] = [motif_symbol.split('.')[0]]

    return(motif_dict)

def addTFOrthoInfo(motifDictionary, orthology_info):
    with open(orthology_info) as TF_orthology_info:
        for line in csv.reader(TF_orthology_info, delimiter = '\t'):
            try:
                motifDictionary[line[14]].append(line[:6] + line[13:])
            except KeyError:
                pass
    return(motifDictionary)

def dupSortOut(motifDictionary):
    duplicate_dict = OrderedDict()
    for motif_id, duplicate_info in motifDictionary.items():
        if len(duplicate_info) == 1:
            pass
            # print(motif_id + '\t' + duplicate_info[0] + '\tNo ortho group info...')
        if len(duplicate_info) == 2:
            pass
            # print(motif_id + '\t' + str(duplicate_info))
        if len(duplicate_info) > 2:
            # print(motif_id + '\t' + duplicate_info[0] + '\tDuplicates present...')
            for duplicate_number, duplicate_orthogroup in enumerate(duplicate_info[1:], start = 1):
                # print('\t' + motif_id + '\t' + str(duplicate_orthogroup))
                duplicate_dict[duplicate_info[0] + '_v' + str(duplicate_number)] = [motif_id, duplicate_orthogroup]

    for duplicated_TF in set([item[0] for item in duplicate_dict.values()]):
        del motifDictionary[duplicated_TF]

    newMotifDictionary = {}
    for gene_symbol, ortho_info in motifDictionary.items():
        newMotifDictionary[ortho_info[0]] = [gene_symbol, ortho_info[1:]]

    for duplicated_motif, ortho_info in duplicate_dict.items():
        newMotifDictionary[duplicated_motif] = ortho_info

    for motif, info in newMotifDictionary.items():
        if len(info[1]) == 0:
            info[1].append(['NA', 'NULL', 'NULL', 'NULL', 'NULL', 'NULL', 'NA', 'NA', 'NA'])

    return(OrderedDict(sorted(newMotifDictionary.items())))


# cock = dupSortOut(addTFOrthoInfo(parseJASPAR(sys.argv[1]), sys.argv[4]))

# print(set([item[0] for item in cock.values()]))


# for ip, bip in cock.items():
#     print(ip + '\t' + str(bip))
#
# sys.exit()

def addEnsToJASPAR(motifDictionary, hsapOrthos, mmusOrthos):
    with open(hsapOrthos) as hsap_ortho_info:
        for line in csv.reader(hsap_ortho_info, delimiter = '\t'):
            try:
                motifDictionary[line[0]].append(line[1:])
            except KeyError:
                motifDictionary[line[0]+'_v1'].append(line[1:])
                motifDictionary[line[0]+'_v2'].append(line[1:])

            try:
                motifDictionary[line[0]+'_v3'].append(line[1:])
            except KeyError:
                pass
            try:
                motifDictionary[line[0]+'_v4'].append(line[1:])
            except KeyError:
                pass

    for motif, info in motifDictionary.items():
        # print(motif + '\t\t' + str(info))
        if len(info) == 2:
            info.append(['NA', 'NA', 'NA'])


    with open(mmusOrthos) as mmus_ortho_info:
        for line in csv.reader(mmus_ortho_info, delimiter = '\t'):
            try:
                motifDictionary[line[0]].append(line[1:])
            except KeyError:
                motifDictionary[line[0]+'_v1'].append(line[1:])
                motifDictionary[line[0]+'_v2'].append(line[1:])

            try:
                motifDictionary[line[0]+'_v3'].append(line[1:])
            except KeyError:
                pass
            try:
                motifDictionary[line[0]+'_v4'].append(line[1:])
            except KeyError:
                pass


        for motif, info in motifDictionary.items():
            if len(info) == 3:
                info.append(['NA', 'NA', 'NA'])

        # print(motif + '\t\t' + str(info))


    return(motifDictionary)

def getMatrixQualityResultsVals(motifDictionary, JASPAR_restultsDir, species):
    whole_set = {}
    for motif, motif_info in tqdm(motifDictionary.items(), desc = 'Reading ' + species + ' matrix_quality results files'):
        whole_set[motif] = 'NA'

        motif_results_dict = OrderedDict()
        try:
            for results_file in sorted(os.listdir(JASPAR_restultsDir + motif_info[0].replace('_V', '_v')), reverse=True):
                # print ('pop')
                if species + 'scan_positive_set_score_distrib.tab' in results_file:

                    with open(JASPAR_restultsDir + motif_info[0].replace('_V', '_v') + '/' +  results_file) as permuted_matix_prob_dist:
                        for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                            if line[0][0] != '#':
                                motif_results_dict[line[0]] = [line[4]]

                if species + 'scan_positive_set_1000perm' in results_file:

                    with open(JASPAR_restultsDir + motif_info[0].replace('_V', '_v') + '/' +  results_file) as permuted_matix_prob_dist:
                        for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                            if line[0][0] != '#':
                                if line[0] in motif_results_dict:
                                    motif_results_dict[line[0]].append(line[4])

                if species + 'matrix.tf' in results_file:
                    matrix_path = os.path.join(JASPAR_restultsDir, motif_info[0].replace('_V', '_v'), results_file)

        except FileNotFoundError:
            try:
                for results_file in sorted(os.listdir(JASPAR_restultsDir + motif_info[0].title()), reverse=True):
                    if species + 'scan_positive_set_score_distrib.tab' in results_file:

                        with open(JASPAR_restultsDir + motif_info[0].title() + '/' +  results_file) as permuted_matix_prob_dist:
                            for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                if line[0][0] != '#':
                                    motif_results_dict[line[0]] = [line[4]]

                    if species + 'scan_positive_set_1000perm' in results_file:

                        with open(JASPAR_restultsDir + motif_info[0].title() + '/' +  results_file) as permuted_matix_prob_dist:
                            for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                if line[0][0] != '#':
                                    if line[0] in motif_results_dict:
                                        motif_results_dict[line[0]].append(line[4])

                    if species + 'matrix.tf' in results_file:
                        matrix_path = os.path.join(JASPAR_restultsDir, motif_info[0].title(), results_file)

            except FileNotFoundError:
                try:
                    for results_file in sorted(os.listdir(JASPAR_restultsDir + string.capwords(motif_info[0])), reverse=True):
                        if species + 'scan_positive_set_score_distrib.tab' in results_file:

                            with open(JASPAR_restultsDir + string.capwords(motif_info[0]) + '/' +  results_file) as permuted_matix_prob_dist:
                                for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                    if line[0][0] != '#':
                                        motif_results_dict[line[0]] = [line[4]]

                        if species + 'scan_positive_set_1000perm' in results_file:

                            with open(JASPAR_restultsDir + string.capwords(motif_info[0]) + '/' +  results_file) as permuted_matix_prob_dist:
                                for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                    if line[0][0] != '#':
                                        if line[0] in motif_results_dict:
                                            motif_results_dict[line[0]].append(line[4])

                        if species + 'matrix.tf' in results_file:
                            matrix_path = os.path.join(JASPAR_restultsDir, string.capwords(motif_info[0]), results_file)

                except FileNotFoundError:
                    try:
                        for results_file in sorted(os.listdir(JASPAR_restultsDir + motif_info[0].lower()), reverse=True):
                            if species + 'scan_positive_set_score_distrib.tab' in results_file:

                                with open(JASPAR_restultsDir + motif_info[0].lower() + '/' +  results_file) as permuted_matix_prob_dist:
                                    for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                        if line[0][0] != '#':
                                            motif_results_dict[line[0]] = [line[4]]

                            if species + 'scan_positive_set_1000perm' in results_file:

                                with open(JASPAR_restultsDir + motif_info[0].lower() + '/' +  results_file) as permuted_matix_prob_dist:
                                    for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                        if line[0][0] != '#':
                                            if line[0] in motif_results_dict:
                                                motif_results_dict[line[0]].append(line[4])

                            if species + 'matrix.tf' in results_file:
                                matrix_path = os.path.join(JASPAR_restultsDir, motif_info[0].lower(), results_file)

                    except FileNotFoundError:
                        for results_file in sorted(os.listdir(JASPAR_restultsDir + string.capwords(motif_info[0].split('::')[0]) + '::' + string.capwords(motif_info[0].split('::')[1])), reverse=True):
                            if species + 'scan_positive_set_score_distrib.tab' in results_file:

                                with open(JASPAR_restultsDir + string.capwords(motif_info[0].split('::')[0]) + '::' + string.capwords(motif_info[0].split('::')[1]) + '/' +  results_file) as permuted_matix_prob_dist:
                                    for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                        if line[0][0] != '#':
                                            motif_results_dict[line[0]] = [line[4]]

                            if species + 'scan_positive_set_1000perm' in results_file:

                                with open(JASPAR_restultsDir + string.capwords(motif_info[0].split('::')[0]) + '::' + string.capwords(motif_info[0].split('::')[1]) + '/' +  results_file) as permuted_matix_prob_dist:
                                    for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                        if line[0][0] != '#':
                                            if line[0] in motif_results_dict:
                                                motif_results_dict[line[0]].append(line[4])

                            if species + 'matrix.tf' in results_file:
                                matrix_path = os.path.join(JASPAR_restultsDir, string.capwords(motif_info[0].split('::')[0]) + '::' + string.capwords(motif_info[0].split('::')[1]), results_file)

        for weight_value, scores in motif_results_dict.items():
            # print(weight_value + str(scores))
            if len(scores) > 1:
                scores.append(float(scores[0]) - float(scores[1]))
                scores.append((scores[2]/float(scores[1])) * 100 )
                # print(inky + '\t' + str(pinky))
                if scores[3] >= 5:
                    if whole_set[motif] == 'NA':
                        whole_set[motif] = (matrix_path, scores[0])
                    else:
                        break
        if whole_set[motif] == 'NA':
            whole_set[motif] = (matrix_path, 'NA')
    # for tf, pval in whole_set.items():
    #     print(tf + '\t' + str(pval))

    return(whole_set)
    # sys.exit()
        # print('\t' + motif + '\tgets\t' + JASPAR_restultsDir + motif_info[0].title())




def getCSPSSMsMatrixQualityResultsVals(CSPSSMdir):
    motif_results_dict = OrderedDict()
    try:
        for results_file in sorted(os.listdir(CSPSSMdir), reverse=True):
            # print(results_file)
            if '__' not in results_file:
                if 'scan_positive_set_score_distrib.tab' in results_file:
                    # print (results_file)

                    with open(os.path.join(CSPSSMdir,  results_file)) as permuted_matix_prob_dist:
                        for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                            if line[0][0] != '#':
                                motif_results_dict[line[0]] = [line[4]]

                if 'scan_positive_set_1000perm' in results_file:
                    # print (results_file)

                    with open(os.path.join(CSPSSMdir, results_file)) as permuted_matix_prob_dist:
                        for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                            if line[0][0] != '#':
                                if line[0] in motif_results_dict:
                                    motif_results_dict[line[0]].append(line[4])

                if 'matrix.tf' in results_file:
                    matrix_path = os.path.join(CSPSSMdir, results_file)

        for weight_value, scores in motif_results_dict.items():
            # print(weight_value + str(scores))
            if len(scores) > 1:
                scores.append(float(scores[0]) - float(scores[1]))
                scores.append((scores[2]/float(scores[1])) * 100 )
                # print(weight_value + str(scores))
                if scores[3] >= 5:
                    pval = (matrix_path, scores[0])
                    break
                else:
                    pval = (matrix_path, 'NA')
        # print(pval)
    except FileNotFoundError:
        # print('No file present. pval = NA')
        pval = 'NA'
    return(pval)

def getCWPSSMsMatrixQualityResultsVals(CWPSSMdir):
    pval_list          = []
    try:
        for species in ['mz', 'pn', 'ab', 'nb', 'on']:
            motif_results_dict = OrderedDict()
            for results_file in sorted(os.listdir(CWPSSMdir), reverse=True):
                # print(results_file)
                if '__' not in results_file:
                    if '_' + species + '_scan_positive_set_score_distrib.tab' in results_file:
                        # print (os.path.join(CWPSSMdir,  results_file))

                        with open(os.path.join(CWPSSMdir,  results_file)) as permuted_matix_prob_dist:
                            for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                if line[0][0] != '#':
                                    motif_results_dict[line[0]] = [line[4]]

                    if '_' +  species + '_scan_positive_set_1000perm' in results_file:
                        # print (os.path.join(CWPSSMdir,  results_file))

                        with open(os.path.join(CWPSSMdir, results_file)) as permuted_matix_prob_dist:
                            for line in csv.reader(permuted_matix_prob_dist, delimiter = '\t'):
                                if line[0][0] != '#':
                                    if line[0] in motif_results_dict:
                                        motif_results_dict[line[0]].append(line[4])

                    if species + '_matrix.tf' in results_file:
                        matrix_path = os.path.join(CWPSSMdir,  results_file)

            switch = 0
            for weight_value, scores in motif_results_dict.items():
                # print(weight_value + str(scores))
                if len(scores) > 1:
                    scores.append(float(scores[0]) - float(scores[1]))
                    scores.append((scores[2]/float(scores[1])) * 100 )
                    # print(weight_value + str(scores))
                    if scores[3] >= 5:
                        pval_list.append((matrix_path, scores[0]))
                        switch = 1
                        break
                    else:
                        # pval_list.append('NA')
                        pass
                # print(pval)
            if switch == 0:
                pval_list.append((matrix_path, 'NA'))

    except FileNotFoundError:
        # print('No file present. pval = NA')
        pval_list.append(['NA', 'NA', 'NA', 'NA', 'NA'])
    # print(pval_list)
    return(pval_list)

# getCWPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/hE_138378/')
# sys.exit()
######################################################
######################################################
slop = addEnsToJASPAR(dupSortOut(addTFOrthoInfo(parseJASPAR(sys.argv[1]), sys.argv[4])), sys.argv[2], sys.argv[3])



pop = getMatrixQualityResultsVals(slop, sys.argv[5], '_ON__')
pap = getMatrixQualityResultsVals(slop, sys.argv[5], '_AB__')
pep = getMatrixQualityResultsVals(slop, sys.argv[5], '_NB__')

pup = getMatrixQualityResultsVals(slop, sys.argv[7], '_PN_')
pip = getMatrixQualityResultsVals(slop, sys.argv[6], '_MZ_')
######################################################
######################################################


# MA0018  ['CREB1', ['OGid2677', 'MA0018.2', 'CREB1', 'mz.gene.s149.5', 'pn.gene.s430.6', 'ab.gene.s568.5', 'nb.gene.s60.75'],
#                   ['OGid2970', 'MA0018.2', 'CREB1', 'mz.gene.s82.26', 'pn.gene.s117.32', 'ab.gene.s34.62', 'nb.gene.s70.4'],
#                   ['both', 'CREB1', 'ENSG00000118260'], ['both', 'CREB1', 'ENSMUSG00000025958']]

with open('mat_qual_pvals.out1', 'w') as out_file:

    cich_motifs_list = []

    for ship, shop in slop.items():

        shop.append(1e-4)
        for ink, dink in pip.items():
            if ink == ship:
                shop.append(dink)
        for ink, dink in pup.items():
            if ink == ship:
                shop.append(dink)
        for ink, dink in pap.items():
            if ink == ship:
                shop.append(dink)
        for ink, dink in pep.items():
            if ink == ship:
                shop.append(dink)
        for ink, dink in pop.items():
                if ink == ship:
                    shop.append(dink)

        if len(shop[1]) == 0:
            shop.append(['NA', 'NA', 'NA', 'NA', 'NA'])


            print(str(ship) + '\t' + str(shop))


        else:
            # print(str(ship) + '\t' + str(shop))
            cspssm_results_list = []
            for lol in shop[1]:
                # print(lol)
                if isinstance(lol, list):
                    for ip in lol[1:7] + lol[8:]:
                        # print(ip[0])
                        if ip[0] != 'E':
                            if ip != 'NA':
                                # print ('RG-cich/CSPSSMs/' +  ip + ' horn')
                                cspssm_results_list.append(getCSPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' +  ip))
                                if ip[0] != 'E':
                                    if ip[0] != 'N':
                                        if ip not in cich_motifs_list:
                                            # print(ip)
                                            cich_motifs_list.append(ip)

                else:
                    if lol.startswith('O'):
                        pass
                    elif lol.startswith('E'):
                        pass
                    elif lol == shop[0]:
                        pass
                    elif lol == 'NULL':
                        cspssm_results_list.append('NA')
                        # print('No ortholog present...')
                    else:
                        # print('HI')
                        # print ('RG-cich/CSPSSMs/' +  lol + ' dorn')
                        cspssm_results_list.append(getCSPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' +  lol))
                        if lol[0] != 'E':
                            if lol[0] != 'N':
                                if lol not in cich_motifs_list:
                                    # print(lol)
                                    cich_motifs_list.append(lol)

            shop.append(cspssm_results_list)


        # print(shop)
        # print(shop[-1])
        # print(len(shop[-1]))
        # if len(shop[-1]) < 5:
        #     for i in range(5-len(shop[-1])):
        #         shop[-1].append('NA')
        # print(shop)
        # print('\n')

        nu_list = [ship]
        for plink in shop:
            if isinstance(plink, list):
                for dink in plink:
                    if isinstance(dink, list):
                        for sink in dink:
                            nu_list.append(str(sink))
                    elif isinstance(dink, tuple):
                        nu_list.append('|'.join(str(i) for i in dink))
                    else:
                        nu_list.append(str(dink))
            elif isinstance(plink, tuple):
                nu_list.append('|'.join(str(i) for i in plink))
            else:
                nu_list.append(str(plink))

        # print(str(ship) + '\t' + str(shop))
        # print('\t'.join(nu_list))
        out_file.write('\t'.join(nu_list) + '\n')
            # sys.exit()

            # '''
            # here needs to go a loop through shop[1] where the above scanning module needs to be employed on the specific dir for each TF on a species by species manner - then where this yeilds no value, the value needs to be dredged from the ENS tree wide maotifs.
            # '''

print('\n\nOutput file closed.\n\n')

# print(sorted(cich_motifs_list))
# print(len(cich_motifs_list))

with open('cich_PSSMs_in_JAS.out1', 'w') as outer:
    for motif in sorted(cich_motifs_list):

        outer.write(motif+'\n')

####################################################################
e_list = []

for entry in sorted(dir_list):
    if entry not in cich_motifs_list:
        # print(entry)
        e_list.append(entry)

# print('\n\n\n' + str(len(e_list)) + '\n\n\n')

# print(len(cich_motifs_list + e_list))
# for dire in os.listdir('RG-cich/CSPSSMs/'):
#     print(dire)
#
#     if os.path.isdir(dire):
#         print(dire)
#
#         if dire not in cich_motifs_list:
#             print('\t' + dire)

ortho_dict = {}

with open(sys.argv[4]) as ortho_info:
    for line in csv.reader(ortho_info, delimiter = '\t'):
        ortho_dict[line[0]] = line[1:]


with open('mat_qual_pvals.out1') as jas_cspssm_dat, open('mat_qual_pvals_noJAS.out1', 'w') as out_file2:
    for count, line in enumerate(csv.reader(jas_cspssm_dat, delimiter = '\t'), start = 1):

        if line[8] != 'NA':
            out_file2.write('motif_' + str(count) + '\t' + line[0] + ';' + line[1] + ';' + line[8] + '\t' + line[2] + '\t' + ';'.join(line[3:8]) + '\t' + line[17] + '\t' + ';'.join(line[18:23]) + '\t' + ';'.join(line[23:]) + '\n')
        else:
            out_file2.write('motif_' + str(count) + '\t' + line[0] + ';' + line[1] + ';' + line[13] + '\t' + line[2] + '\t' + ';'.join(line[3:8]) + '\t' + line[17] + '\t' + ';'.join(line[18:23]) + '\t' + ';'.join(line[23:]) + '\n')
        st_c = count

    ortho_reclaim       = []
    ortho_reclaim_store = OrderedDict()
    for line in e_list:
        for ing, ping in ortho_dict.items():
            if line in ping:
                if line not in ortho_reclaim:
                    ortho_reclaim.append(ing)
                if line not in ortho_reclaim_store:
                    ortho_reclaim_store[ing] = [ping[12], ping[13], ping[:5]]
                else:
                    ortho_reclaim_store[ing].append([ping[12], ping[13], ping[14], ping[:5]])

    # for ert, dert in ortho_reclaim_store.items():
    #     print(ert + '\t' + str(dert))


    for count, (og, line) in enumerate(sorted(ortho_reclaim_store.items(), key = lambda x: int(x[0].split('_')[0].replace('OG', ''))), start = st_c + 1):
        out_file2.write('motif_' + str(count) + '\tNA;' + str(line[1]) + ';' + '\t'.join([line[0], og, ';'.join(line[2]) , '0.0001' , ';'.join(['NA', 'NA', 'NA', 'NA', 'NA']) , ';'.join(['|'.join(str(i) for i in getCSPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' +  x)) if getCSPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' +  x) != 'NA' else getCSPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' +  x) for x in line[2]])]) + '\n')


def addInCWPSSMvals (resultsMatrix):
    nu_dict = OrderedDict()
    with open(resultsMatrix) as results_to_add_to:
        for line in csv.reader(results_to_add_to, delimiter = '\t'):
            result     = getCWPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/' + line[1].split(';')[2].replace('ENSG00000', 'hE_'))
            new_result = []
            if isinstance(result[0], list):
                for item in result[0]:
                    new_result.append(item)
            else:
                for item in result:
                    if item != 'NA':
                        new_result.append('|'.join(str(i) for i in item))

            line.append(';'.join(new_result))

            nu_dict[line[0]] = line[1:]

    return(nu_dict)
            # print(line[-1])

            # print(';'.join(line[-1]))

parf = addInCWPSSMvals('mat_qual_pvals_noJAS.out1')

# for art, dart in parf.items():
#     print(art + '\t' + '\t'.join(dart) + '\n')



with open('mat_qual_pvals_ALL.out2', 'w') as out_file3:
    for art, dart in parf.items():
        out_file3.write(art + '\t' + '\t'.join(dart) + '\n')

def bigFIMOrun (finalMotifDict):

    print('Making Output architecture...')
    try:
        os.mkdir(os.path.join(sys.argv[8], 'bigFIMOtings'))

        for line , pine in finalMotifDict.items():
            try:
                os.mkdir(os.path.join(sys.argv[8], 'bigFIMOtings', line))

            except FileExistsError:
                print('Dir ' + os.path.join(sys.argv[8], 'bigFIMOtings', line) + ' already exists, writing there... ')

            print(os.path.join(sys.argv[8], 'bigFIMOtings', line, line) + '_meta_info.txt')

            print(line)

            with open(os.path.join(sys.argv[8], 'bigFIMOtings', line, line) + '_meta_info.txt', 'w') as summary_output:
                summary_output.write('\tJASPAR motif:\t' + pine[0].split(';')[0] + '\n')
                summary_output.write('\tGene symbol:\t' + pine[0].split(';')[1] + '\n')
                summary_output.write('\tEnsembl id:\t' + pine[0].split(';')[2] + '\n\n')
                summary_output.write('\tOrthoGroup:\t' + pine[1] + '\n')
                summary_output.write('\t\t' + '\t'.join(pine[2].split(';')) + '\n')
                summary_output.write('\tGeneric Pval:\t' + pine[3] + '\n')
                summary_output.write('\tJASPAR Specific Matracies:\n')
                for jas in pine[4].split(';'):
                    summary_output.write('\t\t' + '\t'.join(jas.split('|')) + '\n')
                summary_output.write('\n')
                summary_output.write('\tCichlid Specific Matracies:\n')
                for cic in pine[5].split(';'):
                    summary_output.write('\t\t' + '\t'.join(cic.split('|')) + '\n')
                summary_output.write('\n')
                summary_output.write('\tCichlid Wide Matracies:\n')
                for ciw in pine[6].split(';'):
                    summary_output.write('\t\t' + '\t'.join(ciw.split('|')) + '\n')
                summary_output.write('\n')

            for sp in ['mz', 'pn', 'ab', 'nb', 'on']:
                try:
                    os.mkdir(os.path.join(sys.argv[8], 'bigFIMOtings', line, sp))

                    for sub_folder in ['1_JASPAR', '2_CS', '3_CW']:
                        os.mkdir(os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder))

                except FileExistsError:
                    print('Dir ' + os.path.join(sys.argv[8], 'bigFIMOtings', line, sp) + ' already exists, writing there... ')


    except FileExistsError:
        print('Dir ' + os.path.join(sys.argv[8], 'bigFIMOtings') + ' already exists, writing there... ')




def promPaths (species):
    if species == 'on':
        promSeqPath = '/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta'
        backgroundPath = 'Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta.0orderMarkovBackgrnd'
    else:
        promSeqPath = '/tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/' + string.capwords(species) + '_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta'
        backgroundPath = string.capwords(species) + '_GeneAnnotation_11092017_FINALcorrected.5kb_promoters.stranded.GENEBED_all.fasta.0orderMarkovBackgrnd'
    return((backgroundPath, promSeqPath))

def pbsWrapper (process, *lines):
    with open(process + '_sub.sh', 'w') as temp_sub_script:
        temp_sub_script.write('#!/bin/bash\n#PBS -N fimo\n#PBS -q Prod\n#PBS -l select=1:mem=100GB:ncpus=1\n#PBS -l walltime=500:00:00\n#PBS -rn\n#PBS -l place=excl\n#PBS -o Scr/logs/' + process + '.STDOUT\n#PBS -e Scr/logs/' + process + '.STDERR\n\nexport PATH="$PATH:/tgac/software/testing/meme/4.11.1/x86_64/bin/"\n\n')
        for line in lines:
            temp_sub_script.write(line + '\n\n')

    subprocess.run(['qsub ' + process + '_sub.sh'], check = True, shell = True)
    time.sleep(.5)
    subprocess.run(['rm -f ' + process + '_sub.sh'], check = True, shell = True)

def transMatrix (inputTransfac, outputMEME):
    transCall = 'transfac2meme ' + inputTransfac + ' > ' + outputMEME
    return (transCall)

def fimoScript (background, promoter, output, matrix, pval):
    fimoCall = 'fimo --bgfile ' + background + ' --o ' + output + ' --thresh ' + pval + ' --max-stored-scores 500000 '  + matrix + ' ' + promoter
    return(fimoCall)



def FIMOrun (finalMotifDict):
    count = 0
    for line , pine in finalMotifDict.items():
        print (line)
        for sub_folder in ['1_JASPAR', '2_CS', '3_CW']:
            for sp in ['mz', 'pn', 'ab', 'nb', 'on']:

            # print(os.path.join(sys.argv[8], 'bigFIMOtings', line, sp))

                print(os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder))
                if sub_folder == '1_JASPAR':
                    for jas in pine[4].split(';'):
                        if jas.split('|')[0] != 'NA':
                            # print(jas.split('|'))
                            if '_' + sp.upper() in jas.split('|')[0] :

                                # print('\t\t' + jas.split('|')[0] + '\t' + pine[3])
                                count = count + 2

                                pbsWrapper(line + '_' + sub_folder.split('_')[1] + '_' + sp, transMatrix(jas.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), jas.split('|')[1]))

                                # sys.exit()

                    #             print(transMatrix(jas.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')) + '\n')
                    #
                    #             print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]) + '\n')
                    #
                    # print('\n')
                    # for jas in pine[4].split(';'):
                    #     if jas.split('|')[0] != 'NA':
                    #         # print(jas.split('|'))
                    #         if jas.split('|')[1] != 'NA':
                    #             if '_' + sp.upper() in jas.split('|')[0] :
                    #                 # print('\t\t' + '\t'.join(jas.split('|')))
                    #                 count = count + 1
                    #
                    #                 print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), jas.split('|')[1]) + '\n')

                    print('\n')

                if sub_folder == '2_CS':
                    for cic in pine[5].split(';'):
                        if cic.split('|')[0] != 'NA':
                            if sp + '.' in cic.split('|')[0]:
                                # print('\t\t' + cic.split('|')[0] + '\t' + pine[3])
                                count = count + 2


                                pbsWrapper(line + '_' + sub_folder.split('_')[1] + '_' + sp, transMatrix(cic.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), cic.split('|')[1]))

                    #             print(transMatrix(cic.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')) + '\n')
                    #
                    #             print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]) + '\n')
                    #
                    # print('\n')
                    # for cic in pine[5].split(';'):
                    #     if cic.split('|')[0] != 'NA':
                    #         # print(cic.split('|'))
                    #         if cic.split('|')[1] != 'NA':
                    #             if sp + '.' in cic.split('|')[0]:
                    #                 # print('\t\t' + '\t'.join(cic.split('|')))
                    #                 count = count + 1
                    #
                    #                 print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), cic.split('|')[1]) + '\n')

                    print('\n')

                if sub_folder == '3_CW':
                    for ciw in pine[6].split(';'):
                        if ciw.split('|')[0] != 'NA':
                            if sp + '_' in ciw.split('|')[0]:
                                # print('\t\t' + ciw.split('|')[0] + '\t' + pine[3])
                                count = count + 2


                                pbsWrapper(line + '_' + sub_folder.split('_')[1] + '_' + sp, transMatrix(ciw.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]), fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), ciw.split('|')[1]))

                    #             print(transMatrix(ciw.split('|')[0], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme')) + '\n')
                    #
                    #             print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'default'), os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), pine[3]) + '\n')
                    #
                    # print('\n')
                    # for ciw in pine[6].split(';'):
                    #     if ciw.split('|')[0] != 'NA':
                    #         # print(cic.split('|'))
                    #         if ciw.split('|')[1] != 'NA':
                    #             if sp + '_' in ciw.split('|')[0]:
                    #                 # print('\t\t' + '\t'.join(cic.split('|')))
                    #                 count = count + 1
                    #
                    #                 print(fimoScript(promPaths(sp)[0], promPaths(sp)[1], os.path.join(sys.argv[8], 'bigFIMOtings', line, sp, sub_folder, 'matQualPval'), os.path.join(sys.argv[8], 'bigFIMOtings', line + '_' + sub_folder.split('_')[1] + '_' + sp + '.meme'), ciw.split('|')[1]) + '\n')

                    print('\n')
    print('\n\n\n' + str(count) + '\n\n\n')





        # print('\tJASPAR motif:\t' + pine[0].split(';')[0])
        # print('\tGene symbol:\t' + pine[0].split(';')[1])
        # print('\tEnsembl id:\t' + pine[0].split(';')[2] + '\n')
        # print('\tOrthoGroup:\t' + pine[1])
        # print('\t\t' + '\t'.join(pine[2].split(';')))
        # print('\tGeneric Pval:\t' + pine[3] + '\n')
        # print('\tJASPAR Specific Matracies:')
        # for jas in pine[4].split(';'):
        #     if js == 'NA':
        #         print('no fimo run for you fuc boi')
        #     else:
        #         print('\t\t' + '\t'.join(jas.split('|')))
        #         fimo --bgfile Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta.0orderMarkovBackgrnd --o on-5kbprom_pvalJASPAR_BG --thresh [input_p-val] --max-stored-scores 500000 [matrix_path] /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta
        #
        #
        # print('\tCichlid Specific Matracies:')
        # for cic in pine[5].split(';'):
        #     print('\t\t' + '\t'.join(cic.split('|')))
        #
        # print('\tCichlid Wide Matracies:')
        # for ciw in pine[6].split(';'):
        #     print('\t\t' + '\t'.join(ciw.split('|')))



        # for element in pine:
        #     if isinstance(element.split(';'), list):
        #         for chunk in element.split(';'):
        #             if isinstance (chunk.split(':'), list) or chunk.split('|') == 'NA':
        #                 print('\t\t' + '\t'.join(chunk.split('|')))
        #             else:
        #                 # print('\t\t' + '\t'.join(element.split(';')))
        #                 break
        #     else:
        #         print('\t' + element)
#
# bigFIMOrun(parf)
# FIMOrun(parf)




# subprocess.run(['rm -rf ' + os.path.join(sys.argv[8], 'bigFIMOtings')], check = True, shell = True)






# def submit2FIMO (pvalDict):

    # fimo --bgfile Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta.0orderMarkovBackgrnd --o on-5kbprom_pvalJASPAR_BG --thresh [input_p-val] --max-stored-scores 500000 [matrix_path] /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/1stexonBLAT/Oreochromis_niloticus.BROADON1.longestCDS_new.5kb_promoters.stranded.GENEBED_all.fasta










# getCWPSSMsMatrixQualityResultsVals('RG-cich/CSPSSMs/hE_138378/')







#
# with open('mat_qual_pvals.out') as jas_cspssm_dat:
#     for count, line in enumerate(csv.reader(jas_cspssm_dat, delimiter = '\t'), start = 1):
#         print('motif_' + str(count) + '\t' + line[0] + '\t' + line[1] + '\t' + '\t'.join(line[2:9]) + '\t' + line[10])
#         st_c = count
#
#     ortho_reclaim       = []
#     ortho_reclaim_store = OrderedDict()
#     for line in e_list:
#         for ing, ping in ortho_dict.items():
#             if line in ping:
#                 if line not in ortho_reclaim:
#                     ortho_reclaim.append(ing)
#                 if line not in ortho_reclaim_store:
#                     ortho_reclaim_store[ing] = [ping[12], ping[13], ping[14], ping[:5]]
#                 else:
#                     ortho_reclaim_store[ing].append([ping[12], ping[13], ping[14], ping[:5]])
#
#     # for ert, dert in ortho_reclaim_store.items():
#     #     print(ert + '\t' + str(dert))
#
#
#     for count, (og, line) in enumerate(sorted(ortho_reclaim_store.items(), key = lambda x: int(x[0].split('_')[0].replace('OG', ''))), start = st_c + 1):
#         print('motif_' + str(count) + '\tNA\t' + str(line[1]) + '\t' + '\t'.join([og, '\t'.join(line[3]), line[0], line[2]]) )


































#
# def createMotifDict(target):
#     motif_dict = {}
#     fail_list  = []
#     # with open(target) as motif_dir:
#     for motif_name in tqdm(os.walk(target)):
#         print(motif_name[0])
#         for fil in motif_name[2]:
#             if '1000perm' in fil:
#                 print('\t' + motif_name[0] + '/' + fil)
#                 print('\t' + str(os.path.getsize(motif_name[0] + '/' + fil)))
#
#                 if os.path.getsize(motif_name[0] + '/' + fil) > 26:
#                     if motif_name[0] not in motif_dict:
#                         motif_dict[motif_name[0]] = [(motif_name[0] + '/' + fil, str(os.path.getsize(motif_name[0] + '/' + fil)))]
#                     else:
#                         motif_dict[motif_name[0]].append((motif_name[0] + '/' + fil, str(os.path.getsize(motif_name[0] + '/' + fil))))
#         try:
#             for imp in motif_dict[motif_name[0]]:
#                 print('\t\t'+'\t'.join(imp))
#         except KeyError:
#             fail_list.append(motif_name[0])
#
#             # motif_dict[ortho[0]] = ortho[1]
#     return (motif_dict, fail_list)
#
# prp, lop = createMotifDict(sys.argv[2])
#
# print(len(prp))
# print(len(lop))
