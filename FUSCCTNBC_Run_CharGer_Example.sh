bp1="BP1GeneList.txt"
pp2="PP2GeneList.txt"
out="charged.demo.tsv"
in="demo.vcf"
pvs1="inheritanceGeneList.txt"
ps1="pathogenicVariants.vcf"
pm5="somaticHotspots.hotspot3d.clusters"

cd /home/mading_md2/Biosoft/CharGer

## GATK
charger -f /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_GATK.BamReadcount.MAFfiltered.VEP.vcf -o /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_GATK.BamReadcount.MAFfiltered.VEP.CharGer.tsv -D --inheritanceGeneList inheritanceGeneList.txt -z emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP.vcf --PP2GeneList inheritanceGeneList.txt 

## Pindel insertions
charger -f /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_Pindel.ins_BamReadcount.MAFfiltered.VEP.vcf -o /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_Pindel.ins_BamReadcount.MAFfiltered.VEP.CharGer.tsv -D --inheritanceGeneList inheritanceGeneList.txt -z emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP.vcf --PP2GeneList inheritanceGeneList.txt 

## Pindel deletions
charger -f /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_Pindel.del_BamReadcount.MAFfiltered.VEP.vcf -o /home/mading_md2/Documents/FUSCCTNBC_Germline_Example/TNBC34_Pindel.del_BamReadcount.MAFfiltered.VEP.CharGer.tsv -D --inheritanceGeneList inheritanceGeneList.txt -z emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP.vcf --PP2GeneList inheritanceGeneList.txt 



