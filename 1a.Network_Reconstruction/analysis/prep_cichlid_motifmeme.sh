#!/bin/sh

############################################################################################

# This script will prepare finalised MEME format PWM for each cichlid motifs
# MEME files of all motif for motif scanning according to confidence levels:
# 2a - cichlid-specific matrices;
# 2b - non-species specific matrices;
# 2c - JASPAR matrices

# These files can be used for fimo motif scanning of noncoding regions

############################################################################################

mkdir /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs
cd /tgac/workarea/group-vh/Tarang/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs


# 1. Create MEME format PWM that can be used for FIMO motif scanning

# Derived from Human
mkdir HumanDerived
cd HumanDerived
head -9 /tgac/scratch/nashw/bigFIMOtings/motif_99/ab/1_JASPAR/motif_99_JASPAR_ab.meme > MEME_header

nano create_memepwm.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 28000 # memory pool for all cores
#SBATCH -t 0-00:45 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

# Mz
cat MEME_header > 2c_JASPAR_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/mz/1_JASPAR/motif_*_JASPAR_mz.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_mz.meme

cat MEME_header > 2a_CS_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/mz/2_CS/motif_*_CS_mz.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_mz.meme

cat MEME_header > 2b_CW_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/mz/3_CW/motif_*_CW_mz.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_mz.meme

# Pn
cat MEME_header > 2c_JASPAR_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/pn/1_JASPAR/motif_*_JASPAR_pn.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_pn.meme

cat MEME_header > 2a_CS_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/pn/2_CS/motif_*_CS_pn.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_pn.meme

cat MEME_header > 2b_CW_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/pn/3_CW/motif_*_CW_pn.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_pn.meme

# Ab
cat MEME_header > 2c_JASPAR_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/ab/1_JASPAR/motif_*_JASPAR_ab.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_ab.meme

cat MEME_header > 2a_CS_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/ab/2_CS/motif_*_CS_ab.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_ab.meme

cat MEME_header > 2b_CW_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/ab/3_CW/motif_*_CW_ab.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_ab.meme

# Nb
cat MEME_header > 2c_JASPAR_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/nb/1_JASPAR/motif_*_JASPAR_nb.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_nb.meme

cat MEME_header > 2a_CS_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/nb/2_CS/motif_*_CS_nb.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_nb.meme

cat MEME_header > 2b_CW_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/nb/3_CW/motif_*_CW_nb.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_nb.meme

# On
cat MEME_header > 2c_JASPAR_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/on/1_JASPAR/motif_*_JASPAR_on.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_on.meme

cat MEME_header > 2a_CS_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/on/2_CS/motif_*_CS_on.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_on.meme

cat MEME_header > 2b_CW_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings/motif_*/on/3_CW/motif_*_CW_on.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_on.meme

# run the above
sbatch create_memepwm.sh

# Derived from Mouse
cd ../
mkdir MouseDerived
cd MouseDerived

head -9 /tgac/scratch/nashw/bigFIMOtings_Mmus/motif_99/ab/1_JASPAR/motif_99_JASPAR_ab.meme > MEME_header

nano create_memepwm.sh

#!/bin/bash -e
#SBATCH -p tgac-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 28000 # memory pool for all cores
#SBATCH -t 0-00:45 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

# Mz
cat MEME_header > 2c_JASPAR_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/mz/1_JASPAR/motif_*_JASPAR_mz.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_mz.meme

cat MEME_header > 2a_CS_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/mz/2_CS/motif_*_CS_mz.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_mz.meme

cat MEME_header > 2b_CW_mz.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/mz/3_CW/motif_*_CW_mz.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_mz.meme

# Pn
cat MEME_header > 2c_JASPAR_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/pn/1_JASPAR/motif_*_JASPAR_pn.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_pn.meme

cat MEME_header > 2a_CS_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/pn/2_CS/motif_*_CS_pn.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_pn.meme

cat MEME_header > 2b_CW_pn.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/pn/3_CW/motif_*_CW_pn.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_pn.meme

# Ab
cat MEME_header > 2c_JASPAR_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/ab/1_JASPAR/motif_*_JASPAR_ab.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_ab.meme

cat MEME_header > 2a_CS_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/ab/2_CS/motif_*_CS_ab.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_ab.meme

cat MEME_header > 2b_CW_ab.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/ab/3_CW/motif_*_CW_ab.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_ab.meme

# Nb
cat MEME_header > 2c_JASPAR_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/nb/1_JASPAR/motif_*_JASPAR_nb.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_nb.meme

cat MEME_header > 2a_CS_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/nb/2_CS/motif_*_CS_nb.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_nb.meme

cat MEME_header > 2b_CW_nb.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/nb/3_CW/motif_*_CW_nb.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_nb.meme

# On
cat MEME_header > 2c_JASPAR_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/on/1_JASPAR/motif_*_JASPAR_on.meme)
tail -q -n +10 "${list[@]}" >> 2c_JASPAR_on.meme

cat MEME_header > 2a_CS_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/on/2_CS/motif_*_CS_on.meme)
tail -q -n +10 "${list[@]}" >> 2a_CS_on.meme

cat MEME_header > 2b_CW_on.meme
list=(/tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/on/3_CW/motif_*_CW_on.meme)
tail -q -n +10 "${list[@]}" >> 2b_CW_on.meme

# run the above
sbatch create_memepwm.sh



# 2. Create text files of information on each motif and its meta information including:
# motifID
# genesymbol
# EnsemblID
# orthoGroup
# cichlid_geneIDs
# Optimal p-values for scanning for either JASPAR, cichlid-wide and cichlid-specific instances of that motif

# Example source location: /tgac/scratch/nashw/bigFIMOtings/motif_99/motif_99_meta_info.txt

# Human Derived
awk 'FNR==1{print "----------------------------------------------------------------------------------"}1' /tgac/scratch/nashw/bigFIMOtings/motif_*/motif_*_meta_info.txt | sed $'s/\t//g' | sed 's/Matracies/Matrices/g' > cichlid_allHs_motifs_meta_info.txt

# Mouse Derived
awk 'FNR==1{print "----------------------------------------------------------------------------------"}1' /tgac/scratch/nashw/bigFIMOtings_Mmus/motif_*/motif_*_meta_info.txt | sed $'s/\t//g' | sed 's/Matracies/Matrices/g' > cichlid_allMm_motifs_meta_info.txt
