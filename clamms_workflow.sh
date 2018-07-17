#! /usr/bin/bash

###########################################################
#
# Bash script for clamms workflow 
#
# By MoBiDic - Version 0.0.1
#
###########################################################

###########################################################
# GLOBAL
###########################################################

RED='\033[0;31m'
LIGHTRED='\033[1;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'


VERSION= "0.0.1"
VERBOSITY= 4


###########################################################
# Help 
###########################################################

help(){
  exit 1
}

###########################################################
# Log
###########################################################

# -- Log variables 

ERROR=1
WARNING=2
INFO=3
DEBUG=4

# -- Log functions 

error() { log ${ERROR} "${RED}ERROR${NC} : $1" ; }
warning() { log ${WARNING} "${YELLOW}WARNING${NC} : $1" ; }
info() { log ${INFO} "${BLUE}INFO${NC} : $1" ; }
debug() { log ${DEBUG} "${LIGHTRED}DEBUG${NC} : $1" ; }

# -- Print log 

echoerr() { echo -e "$@" 1>&2 ; }

log() {

    if [ ${VERBOSITY} -ge $1 ]
    then
      echoerr "[`date +'%Y-%m-%d %H:%M:%S'`] - CromWrap version : ${VERSION} - $2"
    fi
  }


###########################################################
# INSTALL 
###########################################################

#Main steps for CNV calling by clamms (https://github.com/rgcgithub/clamms)
# -- Command : ./clamms_workflow.sh install /PATH/TO/Install

install(){
  debug "install : Installation PATH of clamms is : \"$1\""

  info "install : Checking installation directory of clamms..."

  # - Check if /PATH/TO/Install exists else create directory
  if [ ! -d $1 ] 
  then 
    warning "\"$1\" does not exist. \"$1\" was created"
    mkdir $1
  fi
  info "... Argument Checking : OK"

  info "Starting git clone of CLAMMS repo..."
  cd $1
  git clone https://github.com/rgcgithub/clamms.git
  cd clamms
  info "... Done"

  info "Lauching Make ..."
  make
  info "... Done"
}

###########################################################
# MAPABILITY INSTALL - 1T
###########################################################

# -- Command : ./clamms_workflow.sh mapinstall /PATH/TO/INSTALL /PATH/TO/BigWigToWig

mapinstall(){

  debug "mapinstall : Installation PATH of mapability is : \"$1\""
  debug "mapinstall : BigWigToWig PATH is : \"$2\""

  info "Checking mapinstall's arguments ..."

  # - Check if /PATH/TO/INSTALL exists else create directory
  if [ ! -d $1 ]
  then 
    warning "\"$1\" doens't exist. \"$1\" was created"
    mkdir $1
  fi 
  # - Check if /PATH/TO/bigWigToWig exists
  if [ ! -f $2 ]
  then 
    error "\"$2\" is not a correct PATH to BigWigToWig !"
    help
  fi 
  # - Check if /PATH/TO/bigWigToWig is a PATH to the soft
  if ! [[ "$2" =~ bigWigToWig$ ]]
  then 
    error "\"$2\" is not bigWigToWig ! Please chose a correct PATH."
    help 
  fi 
  info "... Arguments checking done"
  info "Installing Mapability ..."
  cd $1
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
  $2 wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wig
  grep -v '^#' wgEncodeCrgMapabilityAlign100mer.wig | sed 's/^chr//g' > ~/PROJECTS/EXOMES/ressources/mappability.bed
  info "... Done !"

}


###########################################################
# CREATE WINDOWS BED - 1T
###########################################################

# -- Command : ./clamms_workflow.sh windowsBed /PATH/TO/annotate_windows.sh INSER_SIZE

#A faire tourner avec "chr", enlever "chr" en sortie du windows.bed
#le github préconise de faire tourner sans le "chr" => prevoir un genome fasta sans chr
#pour l'insert size prevoir un peu plus long que le fragment size ex: 200pb -> 250pb

windowsBed(){

  CLAMS_DIR=$1
  INSER_SIZE=$2
  
  debug "windowsBed : Clams directory is : \"$1\""
  debug "windows Bed : Inser size is : \"$2\""

  info "Checking windowsBed's arguments..."

  # - Check if /PATH/TO/CLAMS exists
  if [ ! -d $1 ]
  then 
    error "\"$1\" is not a correct path to clams directory"
    help 
  fi 
  # - Check if /PATH/TO/annotate_windows.sh exists 
  if [ ! -f "$1/annotate_windows.sh" ]
  then 
    error "\"annotate_windows.sh\" in \"$1\" does not exist !"
    help 
  fi 
  # - Check if /PATH/TO/data/clams_special_regions.bed exists 
  if [ ! -f "$1/data/clams_special_regions.bed" ]
  then 
    error "\"data/clams_special_regions.bed\" in \"$1\" does not exist !"
    help 
  fi 

  # - Check if INSER_SIZE is an integer
  if ! [[ "$2" =~ ^[0-9]+$ ]]
  then 
    error "\"$2\" is not an integer !"
    help 
  fi 

  info "... Argument checking : done !"
  info "Launching annotate_windows.sh ..."

  chmod +x ${CLAMMS_DIR}/annotate_windows.sh
  ${CLAMMS_DIR}/annotate_windows.sh  ~/resources/hg19/S04380110_Regions_nochr.bed  ~/resources/hg19/ucsc.hg19.nochr.fasta  ~/PROJECTS/EXOMES/ressources/mappability.bed ${INSERT_SIZE} ${CLAMMS_DIR}/data/clamms_special_regions.bed > ~/PROJECTS/EXOMES/ressources/windows_nochr_S04380110_${INSERT_SIZE}pb.bed

  info "... annotate_windoxs.sh done !"
}



###########################################################
# NORMALIZE COVERAGE DATA
###########################################################

# -- After Samtools 

normalize(){

  # - Command ./clams.sh normalize /PATH/TO/coverage /PATH/TO/clams_dir SAMPLE

  COVERAGE_PATH=$1
  CLAMS_DIR=$2
  SAMPLE=$3

  debug "normalize : Coverage path is : \"$1\""
  debug "normalize : Clams directory is : \"$2\""
  debug "normalize : Sample name is : \"$3\""
 
  info "Checking normalize's arguments ..."

  # - Check if Coverage Path exists 
  if [ ! -f ${COVERAGE_PATH} ]
  then
    error "\"${COVERAGE_PATH}\" does not exist !"
    help 
  fi 
  
  # - Check 
  if [ ! -d ${CLAMS_DIR} ]
  then 
    error "\"${CLAMS_DIR}\" does not exist ! "
    help 
  fi 

  info "... Argument checking : done !"
  info "Lauching normalize ..."

  cd ${COVERAGE_PATH}
  ls *.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE ; do  ${CLAMMS_DIR}/normalize_coverage ${SAMPLE}.coverage.nochr.bed windows.bed >${SAMPLE}.norm.coverage.bed ;done

  info "... normalize Done !"
}


###########################################
# PICARDMERGED - TO DO 
###########################################

#create metrics matrix for KD tree, all in one: (be careful, column position of picard metrics could change, better to do a grep instead)
metricsMatrix(){

  # - Command ./clams.sh metricsMatrix /PATH/TO/hs_metrics SAMPLE


  HS_FOLDER=$1
  SAMPLEID=$2

  debug "metricsMatrix : HS_FOLDER is : \"$1\""
  debug "metricsMatrix : SAMPLEID is : \"$2\""

  info "Checking metricsMatrix's arguments ..."

  # - Check if HS_FOLDER exists 
  if [ ! -d ${HS_FOLDER} ]
  then 
    error "\"${HS_FOLDER}\" does not exist !"
    help 
  fi 

  info "... Argument checking : done !"
  info "Launching metricsMatrix ..."
  
  cd ${HS_FOLDER}

  for i in *hs_metrics.txt ; do 
    echo "SAMPLE"> ${i/hs_metrics/kd_sample} ;
    echo ${i/_hs_metrics.txt} >> ${i/hs_metrics/kd_sample}; 
    grep -v "#" $i | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${i/hs_metrics/kd_hsmetrics} ;
    grep -v "#" ${i/hs/insertsize} | head -3 | tail -2 | cut -f5 > ${i/hs/kd_insertsize} ;
    paste ${i/hs_metrics/kd_sample}  ${i/hs_metrics/kd_hsmetrics} ${i/hs_/kd_insertsize_} > ${i/hs_/kdTree_}     ; 
  done

  #for 1 sample
  #how to get sample ID  argument?
  echo "SAMPLE" > ${SAMPLEID}_kd_sample.txt
  echo ${SAMPLEID} >> ${SAMPLEID}_kd_sample.txt 
  #path to multiple metrics
  grep -v "#" ${SAMPLEID}_hs_metrics.txt | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${SAMPLEID}_kd_hsmetrics.txt 
  #path to insertsize metrics
  grep -v "#" ${SAMPLEID}_insertsize_metrics.txt | head -3 | tail -2 | cut -f5 > ${SAMPLEID}_kd_insertsize_metrics.txt 
  paste ${SAMPLEID}_kd_sample.txt ${SAMPLEID}_kd_hsmetrics.txt   ${SAMPLEID}_kd_insertsize_metrics.txt >  ${SAMPLEID}_kdTree_metrics.txt 
  # $SAMPLEID_kdTree_metrics.txt a ecrire dans le répertoire kdTreeMetrics

  #head the file with parameter names
  cat ${SAMPLEID}_kdTree_metrics.txt  > ALL_kdTreeMetrics.txt
  #fill with the data (check if conversion of "," into "." is nedded)
  for i in *kdTree_metrics.txt; do t -1 $i | sed 's/,/./g' >> ALL_kdTreeMetrics.txt ; done

  #then move to répertoire kdTreeMetrics 
  cp ${SAMPLEID}_kdTree_metrics.txt kdTreeMetricsDir/

  info "... metricsMatrix done !"

}


###########################################
# MAKEKDTREE
###########################################

makekdtree(){

  # - Command ./clams.sh makekdtree knn outdir 
  

  #############################
  # TO DO 
  #############################
  # ajouter en option le outputDir (à tester) 
  # corriger le Rscript pour ne traiter que le sample

  KNN=$1 
  KD_OUT=$2
  
  debug "makekdtree : KNN is  \"$1\""
  debug "makekdtree : KD_OUT is \"$2\"" 

  info "Checking makekdtree's arguments ..."

  if [ ! -d ${KD_OUT} ]
  then 
    error "\"${KD_OUT}\" does not exist !"
    help
    # Rajouter une ligne pour créer le répertoire s'il n'est pas saisi en dur ? 
  fi 

  if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
  then 
    error "\"${KNN}\" is not a correct value ! Must be an integer."
    help 
  fi 

  info "... Argument checking : done !"

  info "Launching KdTree RScript ..."

  
  Rscript ~/PROJECTS/EXOMES/compute_kdTree.Rscript ${KNN} ALL_kdTreeMetrics.txt ${KD_OUT}

  #sort nns.txt for the next JOIN
  for i in *.${KNN}nns.txt
  do 
    sort $i > ${i/txt/sort.txt}
  done

  info "... KdTree RScript done !"
}


############################################
# MODEL + CNV CALLING 
############################################

cnvCalling(){
  #determine sex  for each sample
  ls *.norm.cov.bed | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.txt
  #include absolute path
  ls *.norm.cov.bed | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$PWD/$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.txt

  #include absolute path and sort
  ls *.norm.cov.bed |sort | while read FILE; do SAMPLE=`echo "$FILE" | cut -d '.' -f 1` ; echo -e -n "$SAMPLE\t$PWD/$FILE\t" ;  grep "Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }' ; done > sample.file.sex.sort.txt

  #sort sample sex file for the next JOIN
  sort sample.file.sex.txt > sample.file.sex.sort.txt

  #same call excluding family with grep with argument $FAMILY_GREP, here for the example B00IX1C and B00IX1B, with nochr
  ls *.norm.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 3` ; join ../../$SAMPLE.${KNN}nns.sort.txt  sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 2- > $SAMPLE.${KNN}nns.ref.panel.sex.txt; grep -vE "$FAMILY_GREP" $SAMPLE.${KNN}nns.ref.panel.sex.txt > $SAMPLE.${KNN}nns.ref.panel.sex.familyExcluded.txt    ; $CLAMMS_DIR/fit_models $SAMPLE.${KNN}nns.ref.panel.sex.familyExcluded.txt ../windows_V6_S07604514_nochr_MERGE.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.nochr.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;    done

  #same call for solo without family grep  with nochr
  ls *.norm.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLE; do  SEX=`echo "$SAMPLE" | join - sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 3` ; join ../../$SAMPLE.${KNN}nns.sort.txt  sample.file.sex.sort.nochr.txt | tr ' ' '\t' | cut -f 2- > $SAMPLE.${KNN}nns.ref.panel.sex.txt;  $CLAMMS_DIR/fit_models $SAMPLE.${KNN}nns.ref.panel.sex.txt ../windows_V6_S07604514_nochr_MERGE.bed >$SAMPLE.models.bed ; $CLAMMS_DIR/call_cnv $SAMPLE.norm.coverage.nochr.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed ;    done


}


############################################
# ANNOTATION 
############################################

annotation(){
  export CI=Sample_name_du_CI
  export DAD=Sample_name_du_pere
  export MUM=Sample_name_du_mere


  #export FAMILY_LIST=`cat $CI_family.text`
  export FAMILY_LIST="$DAD $MUM"   #Argument-liste à donner en debut de run 
  export FAMILY_BED=`sed 's/ /.cnv.bed /g; s/$/.cnv.bed/' <<< $FAMILY_LIST` #prevoir le chemin vers les bed
  export FAMILY_GREP=`sed 's/ /|/g' <<< $FAMILY_LIST`

  #ALL in 2 lines :
  bedtools intersect -a $CI.cnv.bed -b $FAMILY_BED -loj -wao | sort -k1,1 -k2,2n | bedtools merge -c 4,5,6,7,8,9,10,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse -delim "|" > ${CI}.HERIT-DN.bed 

  #add annotated header and annotate
  cat ~/PROJECTS/EXOMES/header.annotated.bed > ${CI}.HERIT-DN.annotated.bed

  #annotate with all datas
  bedtools intersect -a ${CI}.HERIT-DN.bed -b /home/puce/resources/hg19/gene_annot_hg19_final.bed -loj >> ${CI}.HERIT-DN.annotated.bed

  #add size of CNV + link to decipher genome browser
  awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}' ${CI}.HERIT-DN.annotated.bed > ${CI}.HERIT-DN.annotated.final.bed 
}


###########################################################
# MAIN 
###########################################################
