#! /bin/bash


#Main steps for CNV calling by clamms (https://github.com/rgcgithub/clamms)


#INSTALL
git clone https://github.com/rgcgithub/clamms.git
cd clamms
make

export CLAMMS_DIR=~/TOOLS/clamms

############################################
# A faire une seule fois 
############################################
##################################################
##########  STEP 1  ##############################
#create mappability.bed

#get bigWigtowig program
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig

#get mappability bed for 75nt
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
#~/TOOLS/bigWigToWig wgEncodeCrgMapabilityAlign75mer.bigWig wgEncodeCrgMapabilityAlign75mer.wig
#grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' > mappability.bed

#get mappability bed for 100nt
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
~/TOOLS/bigWigToWig wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wig
grep -v '^#' wgEncodeCrgMapabilityAlign100mer.wig | sed 's/^chr//g' > ~/PROJECTS/EXOMES/ressources/mappability.bed



##################################################
##########  STEP 2  ##############################
#create windows bed

#A faire tourner avec "chr", enlever "chr" en sortie du windows.bed

#le github préconise de faire tourner sans le "chr" => prevoir un genome fasta sans chr

#pour l'insert size prevoir un peu plus long que le fragment size ex: 200pb -> 250pb
export INSERT_SIZE=250
chmod +x $CLAMMS_DIR/annotate_windows.sh
$CLAMMS_DIR/annotate_windows.sh  ~/resources/hg19/S04380110_Regions_nochr.bed  ~/resources/hg19/ucsc.hg19.nochr.fasta  ~/PROJECTS/EXOMES/ressources/mappability.bed $INSERT_SIZE $CLAMMS_DIR/data/clamms_special_regions.bed > ~/PROJECTS/EXOMES/ressources/windows_nochr_S04380110_${INSERT_SIZE}pb.bed

#############################################
# Fin une fois 
#############################################



##################################################
##########  STEP 3  ##############################
#compute depths of coverage of bam

#############################################
# FAIT 
#############################################


#Using GATK DepthOfCoverage output files and remove "chr"
#$CLAMMS_DIR/gatk_coverage_to_bed sample.gatk_readDepth_1x_q30.out windows.bed >sample.coverage.bed 
$CLAMMS_DIR/gatk_coverage_to_bed sample.gatk_readDepth_1x_q30.out windows.bed | sed 's/chr//'  >sample.coverage.nochr.bed 



#using samtools
#ls PATH*TO*.bam > bam.list
#using for loop
#for i in bam.list; do samtools bedcov -Q 30 ~/PROJECTS/EXOMES/ressources/windows.bed $i | awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' > ~/PROJECTS/EXOMES/sample_coverage_bed/${i}.coverage.bed; done


#or using paralelization process , modify -P to get number of threads
cat bam.list |xargs -P 1 --max-args 1 $CLAMMS_DIR/compute_coverage.sh



##################################################
##########  STEP 4  ##############################
#normalize coverage data

###########################################
# A faire après le samtoolBedCov
##########################################

#normalize nochr bed
ls *.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE ; do  $CLAMMS_DIR/normalize_coverage $SAMPLE.coverage.nochr.bed windows.bed >$SAMPLE.norm.coverage.bed ;done
#ls *.coverage.bed | cut -d '.' -f 1 | while read SAMPLE ; do  $CLAMMS_DIR/normalize_coverage $SAMPLE.coverage.bed windows.bed >$SAMPLE.norm.cov.bed ;done


#or using parallelized processes
#echo '$CLAMMS_DIR/normalize_coverage $1.coverage.bed windows.bed >$1.norm.cov.bed' \
#>normalize_coverage.sh && chmod +x normalize_coverage.sh
#cat list.of.samples.txt | xargs -P $NUM_PROCESSES --max-args 1 ./normalize_coverage.sh



##################################################
##########  STEP 5  ##############################
#Selecting reference panels using the k-d tre

###########################################
# FAIT
###########################################

#Collect HsMetrics
for i in `cat bam.list` ; do echo $i;java -jar ~/TOOLS/picard.jar CollectHsMetrics I=$i O=${i/.bam/_hs_metrics.txt} R=~/resources/hg19/ucsc.hg19.fasta TARGET_INTERVALS=~/resources/hg19/S04380110_Regions_ucsc_interval.list BAIT_INTERVALS=~/resources/hg19/S04380110_Regions_ucsc_interval.list; done
	
#Collect InsertSizeMetrics
for i in `cat bam.list` ; do echo $i;java -jar ~/TOOLS/picard.jar CollectInsertSizeMetrics I=$i O=${i/.bam/_insertsize_metrics.txt} H=${i/.bam/_histoinsertsize.pdf} M=0.5 ; done


#for i in  *.norm.cov.bed; do echo $i; sed 's/chr//' $i > ${i/bed/nochr.bed} ; done
#for i in *.norm.cov.nochr.bed; do grep -vP "[XY]" $i > ${i/bed/noXY.bed} ; done


#create metrics matrix for KD tree

#for i in *hs_metrics.txt ; do echo "SAMPLE"> ${i/hs_metrics/kd_sample} ;echo ${i/_hs_metrics.txt} >> ${i/hs_metrics/kd_sample}     ; done
#number of columns for hs metrics output seems to have changed#for i in *hs_metrics.txt ; do grep -v "#" $i | head -3 | tail -2 | cut -f50,41,37,21,10,51 > ${i/hs_metrics/kd_hsmetrics} ; done
#for i in *hs_metrics.txt ; do grep -v "#" $i | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${i/hs_metrics/kd_hsmetrics} ; done
#for i in *hs_metrics.txt ; do grep -v "#" ${i/hs/insertsize} | head -3 | tail -2 | cut -f5 > ${i/hs/kd_insertsize} ; done
#for i in *hs_metrics.txt; do paste ${i/hs_metrics/kd_sample}  ${i/hs_metrics/kd_hsmetrics} ${i/hs_/kd_insertsize_} > ${i/hs_/kdTree_} ; done

###########################################
# PICARDMERGED A FAIRE 
###########################################

#create metrics matrix for KD tree, all in one: (be careful, column position of picard metrics could change, better to do a grep instead)
for i in *hs_metrics.txt ; do 
  echo "SAMPLE"> ${i/hs_metrics/kd_sample} ;
  echo ${i/_hs_metrics.txt} >> ${i/hs_metrics/kd_sample}; 
  grep -v "#" $i | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${i/hs_metrics/kd_hsmetrics} ;
  grep -v "#" ${i/hs/insertsize} | head -3 | tail -2 | cut -f5 > ${i/hs/kd_insertsize} ;
  paste ${i/hs_metrics/kd_sample}  ${i/hs_metrics/kd_hsmetrics} ${i/hs_/kd_insertsize_} > ${i/hs_/kdTree_}     ; 
done

#for 1 sample
#how to get sample ID  argument?
echo "SAMPLE" > $SAMPLEID_kd_sample.txt
echo $SAMPLEID >> $SAMPLEID_kd_sample.txt 
#path to multiple metrics
grep -v "#" $SAMPLEID_hs_metrics.txt | head -3 | tail -2 | cut -f51,42,38,21,10,52 > $SAMPLEID_kd_hsmetrics.txt 
#path to insertsize metrics
grep -v "#" $SAMPLEID_insertsize_metrics.txt} | head -3 | tail -2 | cut -f5 > $SAMPLEID_kd_insertsize_metrics.txt 
paste $SAMPLEID_kd_sample.txt $SAMPLEID_kd_hsmetrics.txt   $SAMPLEID_kd_insertsize_metrics.txt >  $SAMPLEID_kdTree_metrics.txt 
# $SAMPLEID_kdTree_metrics.txt a ecrire dans le répertoire kdTreeMetrics

#head the file with parameter names
cat $SAMPLEID_kdTree_metrics.txt  > ALL_kdTreeMetrics.txt
#fill with the data (check if conversion of "," into "." is nedded)
for i in *kdTree_metrics.txt; do t -1 $i | sed 's/,/./g' >> ALL_kdTreeMetrics.txt ; done

#then move to répertoire kdTreeMetrics 
cp $SAMPLEID_kdTree_metrics.txt kdTreeMetricsDir/


#head the file with parameter names
#h -1 `ls *_kdTree_metrics.txt | t -1` > ALL_kdTreeMetrics.txt
#fill with the data (check if conversion of "," into "." is nedded)
#for i in *kdTree_metrics.txt; do t -1 $i | sed 's/,/./g' >> ALL_kdTreeMetrics.txt ; done

#exclude familly members
#grep -v -e "B00IX3B" -e "B00IX3A" ALL_kdTreeMetrics.txt > ALL_kdTreeMetrics_forSAMPLE.txt

###########################################
# FIN PICARDMERGED 
###########################################

#get 32 closest dataset , take ALL_kdTreeMetrics.txt as input (hard coded) => output a SAMPLENAME.32nns.txt file
#use 32 in order to remove family members
#first argument = k-neighbor number to compare  //// second argument = ALL_kdTreeMetrics.txt


###########################################
# A FAIRE : MAKEKDTREE
###########################################


#ajouter en option le outputDir et corriger le Rscript pour ne traiter que le sample

export KNN=50 #Variable à moduler en option 
#Rscript ~/PROJECTS/EXOMES/compute_kdTree.Rscript 50 ALL_kdTreeMetrics.txt
Rscript ~/PROJECTS/EXOMES/compute_kdTree.Rscript $KNN ALL_kdTreeMetrics.txt (+ outputDir)

#sort nns.txt for the next JOIN
for i in *.${KNN}nns.txt; do sort $i > ${i/txt/sort.txt}; done






##################################################
##########  STEP 6  ##############################
#Calling CNVs using the selected reference panels


############################################
# A FAIRE : MODEL + CNV CALLING 
############################################


#determine sex  for each sample
ls *.norm.cov.bed | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.txt
#include absolute path
ls *.norm.cov.bed | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$PWD/$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.txt

#include absolute path and sort
ls *.norm.cov.bed |sort | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$PWD/$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.sort.txt




#sort sample sex file for the next JOIN
sort sample.file.sex.txt > sample.file.sex.sort.txt



#CALL CNVs , include select panel, fit_models and call_cnv
#ls *.norm.coverage.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.txt | tr ' ' '\t' | cut -f 3` ; join $SAMPLE.32nns.txt  sample.file.sex.txt | tr ' ' '\t' | cut -f 2- >$SAMPLE.32nns.ref.panel.sex.txt; $CLAMMS_DIR/fit_models $SAMPLE.32nns.ref.panel.sex.txt windows.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;	done

#same call excluding family with grep , here for the example B00IX1C and B00IX1B
#s *.norm.coverage.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.txt | tr ' ' '\t' | cut -f 3` ; join $SAMPLE.32nns.sort.txt  sample.file.sex.sort.txt | tr ' ' '\t' | cut -f 2- >$SAMPLE.32nns.ref.panel.sex.txt; grep -v -e "B00IX1C" -e "B00IX1B" $SAMPLE.32nns.ref.panel.sex.txt > $SAMPLE.32nns.ref.panel.sex.familyExcluded.txt    ; $CLAMMS_DIR/fit_models $SAMPLE.32nns.ref.panel.sex.familyExcluded.txt windows.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;	done



#with nochr
#ls *.norm.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 3` ; join ../../$SAMPLE.32nns.sort.txt  sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 2- > $SAMPLE.32nns.ref.panel.sex.txt; grep -v -e "B00IX1C" -e "B00IX1B" $SAMPLE.32nns.ref.panel.sex.txt > $SAMPLE.32nns.ref.panel.sex.familyExcluded.txt    ; $CLAMMS_DIR/fit_models $SAMPLE.32nns.ref.panel.sex.familyExcluded.txt ../windows_V6_S07604514_nochr_MERGE.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.nochr.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;    done




#same call excluding family with grep with argument $FAMILY_GREP, here for the example B00IX1C and B00IX1B, with nochr
ls *.norm.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 3` ; join ../../$SAMPLE.${KNN}nns.sort.txt  sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 2- > $SAMPLE.${KNN}nns.ref.panel.sex.txt; grep -vE "$FAMILY_GREP" $SAMPLE.${KNN}nns.ref.panel.sex.txt > $SAMPLE.${KNN}nns.ref.panel.sex.familyExcluded.txt    ; $CLAMMS_DIR/fit_models $SAMPLE.${KNN}nns.ref.panel.sex.familyExcluded.txt ../windows_V6_S07604514_nochr_MERGE.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.nochr.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;    done


#same call for solo without family grep  with nochr
ls *.norm.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 3` ; join ../../$SAMPLE.${KNN}nns.sort.txt  sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 2- > $SAMPLE.${KNN}nns.ref.panel.sex.txt;  $CLAMMS_DIR/fit_models $SAMPLE.${KNN}nns.ref.panel.sex.txt ../windows_V6_S07604514_nochr_MERGE.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.nochr.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;    done



###########################################
# OSEF 
###########################################

################################################
############ OPTIONNAL STEP PCA and PLOT #######

ls *.norm.cov.nochr.bed | while read FILE; do echo -e -n "$FILE\t"; grep "^Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'; done > ref.panel.files.txt 



NUM_SAMPLES=`ls *.norm.cov.bed | wc -l | awk '{print $1}'`
NUM_WINDOWS=`ls *.norm.cov.bed | head -n 1 | xargs awk '$1 != "X" && $1 != "Y" && $NF != 0 {x++} END {print x}'`
#or only with chr10
NUM_WINDOWS=`ls *.norm.cov.bed | head -n 1 | xargs awk '$1 == "chr10" {x++} END {print x}'`
echo -e "$NUM_SAMPLES\t$NUM_WINDOWS" >matrix.txt


#ls *.norm.cov.bed | while read FILE ; do awk '$1 != "X" && $1 != "Y" && $NF != 0 { print $4 }' $FILE | gawk -f $CLAMMS_DIR/transpose.gawk >>matrix.txt ; done
#select only chr10, bug if $NF !=0 is add
ls *.norm.cov.bed | while read FILE ; do awk '$1 == "chr10" { print $4 }' $FILE | gawk -f $CLAMMS_DIR/transpose.gawk >>matrix.txt ; done

~/TOOLS/SVDLIBC/svd -d 4 -o svd-output -r dt matrix.txt
ls *.norm.cov.bed | cut -d '.' -f 1 >sample.names.txt
tail -n +2 svd-output-Ut | tr ' ' '\t' | gawk -f $CLAMMS_DIR/transpose.gawk | paste sample.names.txt - >pca.coordinates.txt

#plotting data with R:
coords <- read.table("pca.coordinates.txt", col.names=c("sample", "pc1", "pc2", "pc3", "pc4"), colClasses=c("character", rep("numeric", 4)))
library(ggplot2)
jpeg("myplot.jpg")
ggplot(coords, aes(x = pc1, y = pc2)) + geom_point() + geom_text(aes(label=sample),hjust=1, vjust=1) + xlim(-0.075, -0.065) + ylim(-0.0115,-0.007)
dev.off()


#select closest norm.cov.bed from 32nns.txt file and concatenate coverage
#add also SAMPLE coverage (not in the 
#create matrix for PCA
for i in `cat SAMPLE.${KNN}nns.txt`; do  j=${i/E498_CP_B00GM/}; echo $j; awk -v pat=$j 'BEGIN{line=""}{line = line"\t"$4}END{print pat,line}' ${i}.norm.cov.nochr.noXY.bed >> matrix_ALL.txt ; done
#for i in *.norm.cov.nochr.noXY.bed; do j=${i/.norm.cov.nochr.noXY.bed/}; j=${j/E498_CP_B00GM/}; echo $j; awk -v pat=$j 'BEGIN{line=""}{line = line"\t"$4}END{print pat,line}' $i >> matrix_ALL.txt ; done





################################################
############ STEP ANNOTATION ###################

############################################
# A FAIRE : ANNOTATION 
############################################


export CI=Sample_name_du_CI
export DAD=Sample_name_du_pere
export MUM=Sample_name_du_mere


#export FAMILY_LIST=`cat $CI_family.text`
export FAMILY_LIST="$DAD $MUM"   #Argument-liste à donner en debut de run 
export FAMILY_BED=`sed 's/ /.cnv.bed /g; s/$/.cnv.bed/' <<< $FAMILY_LIST` #prevoir le chemin vers les bed
export FAMILY_GREP=`sed 's/ /|/g' <<< $FAMILY_LIST`


#ALL in 2 lines :
bedtools intersect -a $CI.cnv.bed -b $FAMILY_BED -loj -wao | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6,7,8,9,10,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse -delim "|" > ${CI}.HERIT-DN.bed 
#bedtools intersect -a $CI.cnv.bed -b $FAMILY_BED -loj -wao -names $FAMILY_LIST | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6,7,8,9,10,19,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse -delim "|" > ${CI}.HERIT-DN.bed 
#bedtools intersect -a $CI.cnv.bed -b $DAD.cnv.bed $MUM.cnv.bed -loj -wao -names pere mere| sort -k1,1 -k2,2n | bedtools merge -c 4,5,6,7,8,9,10,19,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse -delim "|" > ${CI}-test.MP-DN.bed 






#add annotated header and annotate
cat ~/PROJECTS/EXOMES/header.annotated.bed > ${CI}.HERIT-DN.annotated.bed

#annotate with all datas
bedtools intersect -a ${CI}.HERIT-DN.bed -b /home/puce/resources/hg19/gene_annot_hg19_final.bed -loj >> ${CI}.HERIT-DN.annotated.bed


#add size of CNV + link to decipher genome browser
awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}' ${CI}.HERIT-DN.annotated.bed > ${CI}.HERIT-DN.annotated.final.bed 





### Filtre par la qualité des CNV (Q_EXACT >= 0) 
#awk '$10 >= 0 {print $0}' ${CI}.HERIT-DN.annotated.bed > ${CI}.HERIT-DN.annotated.filtered-Q_exact.bed

# Afficher si CNV uniquement homozygote
#awk '$6 == 0 {print $0}' ${CI}-CI.annotated.cnv.sort.filter.bed

#post-traitement du bed original
#cat ~/PROJECTS/EXOMES/header.clamms.bed > $CI.cnv.headered.bed
#cat $CI.cnv.bed >> $CI.cnv.headered.bed



##############################################
#############################################
#           BROUILLON



###### Step 1 optionnel si apparenté disponible
### Ajouter le caractere de novo ou hérité
## Si 2 parents (à tester)
# CNV Intersect avec CI + père + % overlap
bedtools intersect -a $CI.cnv.bed -b $DAD.cnv.bed -wo | awk '{print $0"\tPere_"$NF*100/($3-$2)}' > ${CI}-CI-P.cnv.bed

# CNV Intersect avec CI + mère + % overlap
bedtools intersect -a $CI.cnv.bed -b $MUM.cnv.bed -wo | awk '{print $0"\tMere_"$NF*100/($3-$2)}' > ${CI}-CI-M.cnv.bed

# CNV total herité Mere + père
cat ${CI}-CI-M.cnv.bed ${CI}-CI-P.cnv.bed | sort -k1,1 -k2,2n | bedtools merge -c 38 -o collapse -delim "|" > ${CI}-CI-MP.cnv.bed


#TEST 1 line 
bedtools intersect -a $CI.cnv.bed -b $DAD.cnv.bed $MUM.cnv.bed -loj -wao | sort -k1,1 -k2,2n | bedtools merge -c 24,23,25 -o collapse,collapse,collapse -delim "|" > ${CI}-test.bed


bedtools intersect -a $CI.cnv.bed -b $DAD.cnv.bed $MUM.cnv.bed -loj -wao -names pere mere| sort -k1,1 -k2,2n | bedtools merge -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse,collapse -delim "|" > ${CI}-test.MP-DN.bed 



# CNV de novo + annotation segregation 2 parents

bedtools intersect -a $CI.cnv.bed -b $DAD.cnv.bed -v > ${CI}-CI-ExcludeDADMUM.cnv.bed 
bedtools intersect -a $CI.cnv.bed -b $MUM.cnv.bed -v >> ${CI}-CI-ExcludeDADMUM.cnv.bed 
sort -k1,1 -k2,2n ${CI}-CI-ExcludeDADMUM.cnv.bed | uniq -d | awk '{print $0"\tDeNovo"}'  > ${CI}-DN.cnv.bed

## Si 1 parent (fonctionne)
# CNV de novo + 1 parents 

#bedtools intersect -a $CI.gene.cnv.bed -b $MUM.cnv.bed -wo > $CI-CI-M.cnv.bed
#bedtools intersect -a $CI.gene.cnv.bed -b $CI-CI-M.cnv.bed -v | awk '{print $0"\t""De_novo""\t"$NF*100/($3-$2)}' > $CI-CI-DN.cnv.bed


###### Step 2
### Fusion des fichiers CNV attention si duo ou trio 
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$22"\t"$27"\t"$28"\t"$29"\t"$30"\t"$26}' ${CI}-CI-MP.cnv.bed > ${CI}-CI.cnv.bed
cut -f1,2,3,5,6,7,8,9,10,22,23 ${CI}-CI-DN.cnv.bed >> ${CI}-CI.cnv.bed
awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' ${CI}-CI.cnv.bed

sort -k1,1 -k2,2n  ${CI}-CI.cnv.bed > ${CI}-CI.cnv.sort.bed



cat ${CI}-CI-MP.cnv.bed ${CI}-DN.cnv.bed | sort -k1,1 -k2,2n > ${CI}-CI-MP-DN.cnv.bed





###### Step 3
