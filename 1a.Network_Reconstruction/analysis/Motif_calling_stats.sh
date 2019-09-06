#!/bin/sh

cd /tgac/workarea/Research-Groups/RG-cichlids/matQualOutput

for i in CS_default_*_results.out ; do wc -l $i ; done
for i in CS_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CS_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CS_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CS_matQualPval_*_results.out ; do wc -l $i ; done
for i in CS_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CS_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CS_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CW_default_*_results.out ; do wc -l $i ; done
for i in CW_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CW_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CW_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in CW_matQualPval_*_results.out ; do wc -l $i ; done
for i in CW_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in CW_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in CW_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in JASPAR_default_*_results.out ; do wc -l $i ; done
for i in JASPAR_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in JASPAR_default_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in JASPAR_default_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

for i in JASPAR_matQualPval_*_results.out ; do wc -l $i ; done
for i in JASPAR_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | wc -l ; done
for i in JASPAR_matQualPval_*_results.out ; do cut -f1 $i | sort -u | wc -l ; done
for i in JASPAR_matQualPval_*_results.out ; do awk '($8 < 0.05 )' OFS='\t' $i | cut -f1 | sort -u | wc -l ; done

# get total of TFs (except CW) 
cat CS_default_ab_results.out CS_matQualPval_ab_results.out JASPAR_default_ab_results.out JASPAR_matQualPval_ab_results.out | cut -f1 | sort -u | wc -l
cat CS_default_mz_results.out CS_matQualPval_mz_results.out JASPAR_default_mz_results.out JASPAR_matQualPval_mz_results.out | cut -f1 | sort -u | wc -l
cat CS_default_nb_results.out CS_matQualPval_nb_results.out JASPAR_default_nb_results.out JASPAR_matQualPval_nb_results.out | cut -f1 | sort -u | wc -l
cat CS_default_on_results.out CS_matQualPval_on_results.out JASPAR_default_on_results.out JASPAR_matQualPval_on_results.out | cut -f1 | sort -u | wc -l
cat CS_default_pn_results.out CS_matQualPval_pn_results.out JASPAR_default_pn_results.out JASPAR_matQualPval_pn_results.out | cut -f1 | sort -u | wc -l
# get total of TFs (filtered qval < 0.05)
cat CS_default_ab_results.out CS_matQualPval_ab_results.out JASPAR_default_ab_results.out JASPAR_matQualPval_ab_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_mz_results.out CS_matQualPval_mz_results.out JASPAR_default_mz_results.out JASPAR_matQualPval_mz_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_nb_results.out CS_matQualPval_nb_results.out JASPAR_default_nb_results.out JASPAR_matQualPval_nb_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_on_results.out CS_matQualPval_on_results.out JASPAR_default_on_results.out JASPAR_matQualPval_on_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l
cat CS_default_pn_results.out CS_matQualPval_pn_results.out JASPAR_default_pn_results.out JASPAR_matQualPval_pn_results.out | awk '($8 < 0.05 )' OFS='\t' | cut -f1 | sort -u | wc -l


# GTRD
cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman
wc -l hg38_TFTGTFBS_sites_nr.in
cut -f1 hg38_TFTGTFBS_sites_nr.in | sort -u | wc -l

cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDdata
wc -l mm10_TFTGTFBS_sites_nr.in
cut -f1 mm10_TFTGTFBS_sites_nr.in | sort -u | wc -l


## Looking at numbers after extrapolation
cd /tgac/workarea/Research-Groups/RG-cichlids/GTRDhuman/human_sites/
for i in hg38_tftgtfbs_*07.extrap.output ; do wc -l $i ; done
for i in hg38_tftgtfbs_*07.extrap.output ; do awk -F';' '{print $1}' $i | sort -u | wc -l ; done

cd /tgac/workarea/Research-Groups/RG-cichlids/GTRD_TFTGcoexpCheck
for i in mm10_tftgtfbs_*_GeneAnnotation_11092017_07.extrap.output ; do wc -l $i ; done
wc -l mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.output
for i in mm10_tftgtfbs_*_GeneAnnotation_11092017_07.extrap.output ; do awk -F':' '{print $1}' $i | sort -u | wc -l ; done
awk -F':' '{print $1}' mm10_tftgtfbs_Oreochromis_niloticus_07.extrap.output | sort -u | wc -l
