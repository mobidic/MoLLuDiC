#!/usr/bin/bash

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


VERSION="0.0.1"
VERBOSITY=4


###########################################################
# Help 
###########################################################

help() {
  echo "Voici l'aide"
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
    echoerr "[`date +'%Y-%m-%d %H:%M:%S'`] - MoLLuDiC version : ${VERSION} - $2"
  fi
}


###########################################################
# INSTALL 
###########################################################

#Main steps for CNV calling by clamms (https://github.com/rgcgithub/clamms)
# -- Command : ./clamms_workflow.sh install /PATH/TO/Install

install() {
  
  CLAMMS_DIR=$1
  
  debug "install : Installation PATH of clamms is : \"${CLAMMS_DIR}\""

  info "install : Checking installation directory of clamms..."

  # - Check if /PATH/TO/Install exists else create directory
  if [ ! -d ${CLAMMS_DIR} ] 
  then 
    warning "\"${CLAMMS_DIR}\" does not exist. \"${CLAMMS_DIR}\" was created"
    mkdir ${CLAMMS_DIR}
    help
  fi
  info "... Argument Checking : OK"

  info "Starting git clone of CLAMMS repo..."
  cd ${CLAMMS_DIR}
  git clone https://github.com/rgcgithub/clamms.git
  mv clamms/* .
  rm -rf clamms
  info "... Done"

  info "Lauching Make ..."
  make
  info "... Done"
}

###########################################################
# MAPABILITY INSTALL - 1T
###########################################################

# -- Command : ./clamms_workflow.sh mapinstall /PATH/TO/CLAMMS_DIR /PATH/TO/BigWigToWig

mapinstall() {
  
  CLAMMS_DIR=$1
  INSTALLATION_PATH=${CLAMMS_DIR}/lib4Clamms/hg19
  BW2W_PATH=$2
  
  debug "mapinstall : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "mapinstall : INSTALLATION_PATH of mapability is : \"${INSTALLATION_PATH}\""
  debug "mapinstall : BW2W_PATH is : \"${BW2W_PATH}\""

  info "Checking mapinstall's arguments ..."
  # - Check if /PATH/TO/INSTALL exists else create directory
  if [ ! -d ${INSTALLATION_PATH} ]
  then 
    warning "\"${INSTALLATION_PATH}\" doens't exist. \"${INSTALLATION_PATH}\" was created"
    mkdir ${CLAMMS_DIR}
    mkdir ${CLAMMS_DIR}/lib4Clamms
    mkdir ${CLAMMS_DIR}/lib4Clamms/hg19
  fi 
  # - Check if /PATH/TO/bigWigToWig exists
  if [ ! -f ${BW2W_PATH} ]
  then 
    error "\"${BW2W_PATH}\" is not a correct PATH to BigWigToWig !"
    help
  fi 
  # - Check if /PATH/TO/bigWigToWig is a PATH to the soft
  if ! [[ "${BW2W_PATH}" =~ bigWigToWig$ ]]
  then 
    error "\"${BW2W_PATH}\" is not bigWigToWig ! Please chose a correct PATH."
    help 
  fi 
  info "... Arguments checking done"

  info "Installing Mapability ..."
  cd ${INSTALLATION_PATH}
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
  ${BW2W_PATH} wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wig
  grep -v '^#' wgEncodeCrgMapabilityAlign100mer.wig | sed 's/^chr//g' > mappability.bed
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

  CLAMMS_DIR=$1
  INSERT_SIZE=$2
  INTERVALBEDFILE=$3
  REFFASTA=$4
  #REF = CLAMMS_DIR/lib4CLAMMS/hg19
  CLAMMS_SPECIAL_REGIONS=$5
  #CLAMMS_DIR/data
  LIBRARY=$6
  
  debug "windowsBed : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "windowsBed : INSERT_SIZE is : \"${INSERT_SIZE}\""
  debug "windowsBed : INTERVALBEDFILE is : \"${INTERVALBEDFILE}\""
  debug "windowsBed : REFFASTA is : \"${REFFASTA}\""
  debug "windowsBed : CLAMMS_SPECIAL_REGIONS is : \"${CLAMMS_SPECIAL_REGIONS}\""
  debug "windowsBed : LIBRARY is : \"${LIBRARY}\""

  info "Checking windowsBed's arguments..."
  # - If exists
  ## - Check if /PATH/TO/CLAMS exists
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR}\" is not a correct path to clams directory"
    help 
  fi 
  ## - Check if /PATH/TO/annotate_windows.sh exists 
  if [ ! -f "${CLAMMS_DIR}/annotate_windows.sh" ]
  then 
    error "\"annotate_windows.sh\" in \"${CLAMMS_DIR}\" does not exist !"
    help 
  fi 
  ## - Check if INSERT_SIZE >= 100s 
  if [ ! ${INSERT_SIZE} -ge 100 ]
  then 
    error "\"${INSERT_SIZE}\" have to be greater than 100 !"
    help 
  fi 
  ## - Check if INSER_SIZE is an integer
  if ! [[ "${INSERT_SIZE}" =~ ^[0-9]+$ ]]
  then 
    error "\"${INSERT_SIZE}\" is not an integer !"
    help 
  fi
  ## - Check if INTERVALBEDFILE exists
  if [ ! -f ${INTERVALBEDFILE} ]
  then 
    error "\"${INTERVALBEDFILE}\" does not exist !"
    help 
  fi 
  ## - Check if REFFASTA exists
  if [ ! -f ${REFFASTA} ]
  then 
    error "\"${REFFASTA}\" does not exist !"
    help 
  fi 
  ## - Check if CLAMMS_SPECIAL_REGIONS exists
  if [ ! -f ${CLAMMS_SPECIAL_REGIONS} ]
  then 
    error "\"${CLAMMS_SPECIAL_REGIONS}\" does not exist !"
    help
  fi 
  ## - Check if LIBRARY exists 
  if [ ! -d ${LIBRARY} ]
  then 
    error "\"${LIBRARY}\" does not exist !"
    help 
  fi

  # - Extension checking
  
  if [[ ${INTERVALBEDFILE##\.} != "bed" ]]
  then 
    error "\"${INTERVALBEDFILE}\" must be .bed file !"
  fi 

  if [[ ${REFFASTA##\.} != "fa" ]]
  then
    error "\"${REFFASTA}\" must be .fa file !"
  fi 

  if [[ ${CLAMMS_SPECIAL_REGIONS##\.} != "bed" ]]
  then 
    error "\"${CLAMMS_SPECIAL_REGIONS}\" must be .bed file !"
  fi 
  info "... Argument checking : done !"

  info "Launching annotate_windows.sh ..."
  chmod +x ${CLAMMS_DIR}/annotate_windows.sh
  # - Mkdir if insertSize${INSERT_SIZE} does not exist
  if [ -d ${LIBRARY}/windowsBeds/insertSize${INSERT_SIZE} ]
  then 
    echo "\"${LIBRARY}/windowsBeds/insertSize${INSERT_SIZE}\" does not exist but was created !  "
    mkdir ${LIBRARY}/windowsBeds/insertSize${INSERT_SIZE}
  fi 

  # Ajouter libraire en argument - Check variable MobiDL (_nochr.bed)
  ${CLAMMS_DIR}/annotate_windows.sh ${INTERVALBEDFILE} ${REFFASTA}  ${CLAMMS_DIR}/lib4Clamms/hg19/mappability.bed ${INSERT_SIZE} ${CLAMMS_SPECIAL_REGIONS} > ${LIBRARY}/windowsBeds/insertSize${INSERT_SIZE}/windows_nochr_${INSERT_SIZE}pb.bed

  info "... annotate_windoxs.sh done !"
}



###########################################################
# NORMALIZE COVERAGE DATA
###########################################################

# -- After Samtools 

normalizeFS() {

  # - Command ./clams.sh normalize /PATH/TO/coverage /PATH/TO/clams_dir SAMPLE

  COVERAGE_PATH=$1
  CLAMMS_DIR=$2
  SAMPLEID=$3

  debug "normalize : COVERAGE_PATH is : \"${COVERAGE_PATH}\""
  debug "normalize : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "normalize : SAMPLEID name is : \"${SAMPLEID}\""
 
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
  ls *.coverage.nochr.bed | cut -d '.' -f 1 | while read SAMPLEID ; do  ${CLAMMS_DIR}/normalize_coverage ${SAMPLEID}.coverage.nochr.bed windows.bed > ${SAMPLEID}.norm.coverage.bed ;done

  info "... normalize Done !"
}

normalize(){

  COVERAGE_PATH=$1
  CLAMMS_DIR=$2
  SAMPLEID=$3
  CLAMMSCOVERAGEFILE=$4
  WINDOWS_BED=$5

  debug "normalize : COVERAGE_PATH is : \"${COVERAGE_PATH}\""
  debug "normalize : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "normalize : SAMPLEID name is : \"${SAMPLEID}\""
 
  info "Checking normalize's arguments ..."
  # - If exists
  ## - Check if Coverage Path exists 
  if [ ! -f ${COVERAGE_PATH} ]
  then
    error "\"${COVERAGE_PATH}\" does not exist !"
    help 
  fi
  ## - Check 
  if [ ! -d ${CLAMS_DIR} ]
  then 
    error "\"${CLAMS_DIR}\" does not exist ! "
    help 
  fi 
  
  # - Extension checking 
  if [[ ${CLAMMSCOVERAGEFILE##\.} != "bed" ]]
  then 
    error "\"${CLAMMSCOVERAGEFILE}\" must be .bed file !"
  fi 
 
  if [[ ${WINDOWS_BED##\.} != "bed" ]]
  then
    error "\"${WINDOWS_BED}\" must be .bed file !"
  info "... Argument checking : done !"
  fi

  info "Lauching normalize ..."
  cd ${COVERAGE_PATH}
  ${CLAMMS_DIR}/normalize_coverage ${CLAMMSCOVERAGEFILE} ${WINDOWS_BED} > ${SAMPLEID}.norm.coverage.bed

  info "... normalize Done !"


}


###########################################
# PICARDMERGED - TO DO 
###########################################

#create metrics matrix for KD tree, all in one: (be careful, column position of picard metrics could change, better to do a grep instead)
metricsMatrixFS(){

  # - Command ./clams.sh metricsMatrix /PATH/TO/hs_metrics SAMPLE


  HS_FOLDER=$1

  debug "metricsMatrixFS : HS_FOLDER is : \"$1\""

  info "Checking metricsMatrixFS' arguments ..."

  # - Check if HS_FOLDER exists 
  if [ ! -d ${HS_FOLDER} ]
  then 
    error "\"${HS_FOLDER}\" does not exist !"
    help 
  fi 

  info "... Argument checking : done !"
  info "Launching metricsMatrixFS ..."
  
  cd ${HS_FOLDER}
  ## TOUTES LES VARIABLES SONT A RECUP DANS MOBIDL
  for SAMPLEID in *hs_metrics.txt
  do 
    echo "SAMPLE" > ${SAMPLEID}_kd_sample.txt
    echo ${SAMPLEID} >> ${SAMPLEID}_kd_sample.txt 
    #path to multiple metrics
    grep -v "#" ${SAMPLEID}_hs_metrics.txt | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${SAMPLEID}_kd_hsmetrics.txt 
    #path to insertsize metrics
    # Sample ID insertsizer metrics .txt | sample id hs metrics .txt 
    grep -v "#" ${SAMPLEID}_insertsize_metrics.txt | head -3 | tail -2 | cut -f5 > ${SAMPLEID}_kd_insertsize_metrics.txt 
    paste ${SAMPLEID}_kd_sample.txt ${SAMPLEID}_kd_hsmetrics.txt   ${SAMPLEID}_kd_insertsize_metrics.txt >  ${SAMPLEID}_kdTree_metrics.txt 
    # $SAMPLEID_kdTree_metrics.txt a ecrire dans le répertoire kdTreeMetrics
  done

  #head the file with parameter names
  #cat ${SAMPLEID}_kdTree_metrics.txt  > ALL_kdTreeMetrics.txt
  #fill with the data (check if conversion of "," into "." is nedded)
  for i in *kdTree_metrics.txt; do t -1 $i | sed 's/,/./g' >> ALL_kdTreeMetrics.txt ; done

  #then move to répertoire kdTreeMetrics 
  cp ${SAMPLEID}_kdTree_metrics.txt kdTreeMetricsDir/

  info "... metricsMatrixFS done !"


  #for i in *hs_metrics.txt ; do 
  #  echo "SAMPLE"> ${i/hs_metrics/kd_sample} ;
  #  echo ${i/_hs_metrics.txt} >> ${i/hs_metrics/kd_sample}; 
  #  grep -v "#" $i | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${i/hs_metrics/kd_hsmetrics} ;
  #  grep -v "#" ${i/hs/insertsize} | head -3 | tail -2 | cut -f5 > ${i/hs/kd_insertsize} ;
  #  paste ${i/hs_metrics/kd_sample}  ${i/hs_metrics/kd_hsmetrics} ${i/hs_/kd_insertsize_} > ${i/hs_/kdTree_}     ; 
  #done
}

metricsMatrix() {
  
  HS_FOLDER=$1
  SAMPLEID=$2
  HSMETRICSTXT=$3
  INSERT_SIZE_METRICS_TXT=$4

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
  ## TOUTES LES VARIABLES SONT A RECUP DANS MOBIDL
  echo "SAMPLE" > ${SAMPLEID}_kd_sample.txt
  echo ${SAMPLEID} >> ${SAMPLEID}_kd_sample.txt 
  #path to multiple metrics
  grep -v "#" ${HSMETRICSTXT} | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${SAMPLEID}_kd_hsmetrics.txt 
  #path to insertsize metrics
  # Sample ID insertsizer metrics .txt | sample id hs metrics .txt 
  grep -v "#" ${INSERT_SIZE_METRICS_TXT} | head -3 | tail -2 | cut -f5 > ${SAMPLEID}_kd_insertsize_metrics.txt 
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

if [ $1 == "install" ]
then 
  install $2
fi

if [ $1 == "mapinstall" ]
then 
  mapinstall $2 $3
fi 

if [ $1 == "windowsBed" ]
then 
  windowsBed $2 $3 $4 $5 $6 $7
fi

if [ $1 == "normalizeFS" ]
then 
  normalizeFS $2 $3 $4
fi
