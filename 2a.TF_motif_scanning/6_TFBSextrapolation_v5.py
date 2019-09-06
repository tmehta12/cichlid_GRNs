#!/usr/bin/python
import sys
if len(sys.argv) < 4:
	print('\nNot enough arguments entered.')
	print('Usage: Orthology file (; sep), Promoter FASTA, TFBS set\n')
	sys.exit()

print('\n===================================================================================\n\n')
print('TFBSexptrapolation_v4.py\n\n')

print('\t\t!!! RUN IN PYTHON2.7 !!!\n')
print('!!! Requires pyfasta and biopython and numpy-1.7.0!!!\n\n')

print('Importing...\n')
import json
import time
from pyfasta import Fasta
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

print('Witing to output flie:\t' + sys.argv[3] + '_' + str(sys.argv[2].split('/')[-1].split('.')[0]) + '_' + '07.extrap.output\n')
print('Only excepting matches >=70% identity\n\n')
print('===================================================================================\n\n')

print('Working with orthology file:\t' + sys.argv[1] +'\n')

ref_TG              = []
ref_TG_to_target_TG = {}
target_TG           = {}
with open(sys.argv[1]) as f: #this is the orthology file - here 3 dicts of the ortho info is made

	f.readline()

	for line in f:
		cell             = line.strip().split(";")
		gene_id_ref      = cell[0].split(",")
		gene_id_target   = cell[1].split(",")
		gene_name_ref    = cell[2].split(",")
		gene_name_target = cell[3].split(",")

		for refSp_gene_ensembl_id in gene_id_ref:                             #creates a dict of ens_id: [gene_name] for ref sp
			if refSp_gene_ensembl_id not in ref_TG:
				ref_TG.append(refSp_gene_ensembl_id)

			# for gene_name in gene_name_ref:
			# 	if gene_name.upper() not in ref_TG[ refSp_gene_ensembl_id ]:
			# 		ref_TG[ refSp_gene_ensembl_id ].append(gene_name.upper())
			# 		# print('pap: ' + refSp_gene_ensembl_id +' '+ str(gene_name.upper()))

		# for tarSp_gene_name in gene_name_target:                           #creates a dict of gene_name: [ens_id1, ens_id2...] for target sp
		#
		# 	if tarSp_gene_name.upper() not in target_TG:
		# 		target_TG[ tarSp_gene_name.upper() ] = []
		#
		# 	for tarSP_gene_id in gene_id_target:
		#
		# 		if tarSP_gene_id not in target_TG[ tarSp_gene_name.upper() ]:
		# 			target_TG[ tarSp_gene_name.upper() ].append(tarSP_gene_id)

			# print('\n'+str(target_TG)+'\n')
		for refSp_gene_ensembl_id in gene_id_ref:                               #creates a dict of ens_id_ref: [ens_id_traget_1, ens_id_traget_2...]

			if refSp_gene_ensembl_id not in ref_TG_to_target_TG:
				ref_TG_to_target_TG[ refSp_gene_ensembl_id ] = []

			for tarSP_gene_id in gene_id_target:
				if tarSP_gene_id not in ref_TG_to_target_TG[ refSp_gene_ensembl_id ]:
					ref_TG_to_target_TG[ refSp_gene_ensembl_id ].append(tarSP_gene_id)

# TU_to_gene = {}
# with open("TUSet.txt") as f:
# 	for line in f:
# 		if line[0] == "#":
# 			continue
# 		line = line.strip()
# 		cell = line.split("\t")
# 		TU_to_gene[cell[1]] = []
# 		for n in cell[3].split(","):
# 			n = n.lower()
# 			if n in ref_TG:
# 				for u in ref_TG[n]:
# 					if u not in TU_to_gene[cell[1]]:
# 						TU_to_gene[cell[1]].append(u)
# Uniprot_to_RefSeq = {}
# with open("uniprot2refseq.tab") as f:
# 	f.readline()
# 	for line in f:
# 		line = line.strip()
# 		cell = line.split("\t")
# 		if cell[0] not in Uniprot_to_RefSeq:
# 			Uniprot_to_RefSeq[cell[0]] = []
# 		if len(cell) > 1:
# 			for rf in cell[1].split(";"):
# 				if rf != "":
# 					rf = rf.split(".")[0]
# 					if rf not in Uniprot_to_RefSeq[cell[0]]:
# 						Uniprot_to_RefSeq[cell[0]].append(rf)

print('Loaded.\n')
#MAPPING = json.load(open("MAPPIG.json"))
#PWMs = json.load(open("pwm.json"))
promoter_seq_fasta = Fasta(sys.argv[2]) ######################## Promoter fasta seqences are loaded here!

print('Retreiving TFBS info...')
TFBS = []
with open(sys.argv[3]) as f:

	for line in f:
		# print (line)
		# sys.exit()
		if line[0] == "#":
			continue
		TF     = []
		TARGET = []
		cell   = line.strip().split("\t")

		# print(str(cell[1].split('.')[0]))
		#tf
		if cell[0] not in ref_TG:
			continue

		# for tf_u in ref_TG[cell[0]]:     #this is where the orthology info is searched via TF identity

		if cell[0] in ref_TG_to_target_TG:

			for gene_in_target_sp in ref_TG_to_target_TG[cell[0]]:
				# print('\t' + gene_in_target_sp)
				# time.sleep(.5)
				if gene_in_target_sp not in TF:
					TF.append(gene_in_target_sp)
		# print('\t'.join(TF))
		# time.sleep(1)
# 		#target
		# if cell[1] not in TU_to_gene:          #not needed as above TU_to_gene segment commented out
			# continue
		# for p in TU_to_gene[cell[7]]:
		if cell[1].split('.')[0] in ref_TG_to_target_TG:       #
			for u in ref_TG_to_target_TG[cell[1].split('.')[0]]:
				if u not in TARGET:
					TARGET.append(u)
			# print('\t' + '\t'.join(TARGET))
		else:
			continue
		if len(TF) == 0 or len(TARGET) == 0:
			continue
		####
# 		if len(cell) != 13: # this will change for the cichlid instance !!!!! CHANGE (w/ in_f made)
# 			continue

		BindingSite = ""
		for base in cell[2]:
			if base.isupper():
				BindingSite += base
		# print('\t\t' + BindingSite)
# 		####


		for T in TARGET:
# 			for RefSeq in Uniprot_to_RefSeq[T]: # This change to cich gen names, also below
			if T in promoter_seq_fasta:
				if len(promoter_seq_fasta[T]) > len(BindingSite):
					for a in pairwise2.align.localms(str(promoter_seq_fasta[T]).lower(), BindingSite.lower(), 2, -2, -100, -10):



						if (a[4] - a[3]) >= 7*len(BindingSite)/10:
							 if (a[2]/(a[4] - a[3])) >= 0.7:
								print(format_alignment(*a))
								print('\n')
								print(str(a[4] - a[3]) + ' >= ' + str(7*len(BindingSite)/10))
								print(str(str(a[2]/(a[4] - a[3])) + ' >= 0.7'))
								 # print(str((a[4] - a[3]) > 7*len(BindingSite)/10 ))
								 # print(str((a[2]/(a[4] - a[3])) >= 0.7))
								print('\n')
								HIT           = {}
								HIT['TF']     = TF
								HIT['TARGET'] = T
								# HIT['tfbs']   = a[0][a[3]:a[4]].replace("-", "")
								siteHit = []
								for i in range(len(a[0])):
									if a[1][i] != '-':
										siteHit.append(a[0][i])

								# print(a[1].replace('-',''))
								# print(''.join(siteHit))

								time.sleep(.2)

								HIT['tfbs'] = ''.join(siteHit)
								TFBS.append(HIT)
								# print('\n')
								# for pap, chap in HIT.items():
								# 	print(pap + '\t' + str(chap))
								# print('\n')

print('Done.\n')

with open(sys.argv[3] + '_' + str(sys.argv[2].split('/')[-1].split('.')[0]) + '_' + '07.extrap.output', 'w') as f:
	for t in TFBS:
		for tt in t["TF"]:
			f.write("%s;%s;%s\n" % (tt, t["TARGET"], t["tfbs"]) )

f.close()
print('\nOut file closed.\n')
sys.exit()
