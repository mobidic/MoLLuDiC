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
# CLAMMS DIRECTORIES PREPARATION
###########################################################

dirpreparation() {
  
  OPTION=$1

  debug "dirpreparation : OPTION is : \"${OPTION}\""

  case ${OPTION} in 
    clamms)
      CLAMMS_DIR=$2
      debug "dirpreparation clamms : CLAMMS_DIR is : \"${CLAMMS_DIR}\"" 
      mkdir ${CLAMMS_DIR}
      mkdir ${CLAMMS_DIR}lib4Clamms
      mkdir ${CLAMMS_DIR}lib4Clamms/hg19
      ;;
    library)
      CLAMMS_DIR=$2
      LIBRARY_NAME=$3
      debug "dirpreparation library : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
      debug "dirpreparation library : LIBRARY_NAME is : \"${LIBRARY_NAME}\""
      if [ ! -d ${CLAMMS_DIR}lib4Clamms/hg19 ] 
      then
        error "Invalid Clamms directory format : \"${CLAMMS_DIR}lib4Clamms/hg19\" is missing ! Please use : \"dirpreparation clamms DIRECTORY\" to create correct clamms directories."
        exit 1
      fi
      
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/windowsBeds
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/projects/
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/projects/all
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/projects/all/kdTreeMetrics
      mkdir ${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/projects/all/normCoverageNoChrBeds
      ;;
    size)
      LIBRARY_DIR=$2
      INSERT_SIZE=$3
      debug "dirpreparation size : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
      debug "dirpreparation size : INSERT_SIZE is : \"${INSERT_SIZE}\""
      if  [ ! -d ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds ]
      then 
        error "Invalid Library directory format : \"${LIBRARY_DIR}projects/all/normCoverageNoChrBeds\" is missing ! Please use : \"dirpreparation library CLAMMS_DIR LIBRARY_NAME\" to create correct library directory."
      fi

      if [ -d ${LIBRARY_DIR}windowsBeds ]
      then 
        error "Invalid Library directory format : \"${LIBRARY_DIR}windowsBeds\" is missing ! Please use : \"dirpreparation library CLAMMS_DIR LIBRARY_NAME\" to create correct library directory."
      fi 

      if ! [[ "${INSERT_SIZE}" =~ ^[0-9]+$ ]]
      then 
        error "\"${INSERT_SIZE}\" must be an integer ! Please select a correct value."
        exit 1
      fi 

      mkdir ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}/
      ;;
  esac
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
    mk
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
  INSTALLATION_PATH=${CLAMMS_DIR}lib4Clamms/hg19
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
    mkdir ${CLAMMS_DIR}lib4Clamms
    mkdir ${CLAMMS_DIR}lib4Clamms/hg19
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
# -- Command : ./clamms_workflow.sh windowsBed CLAMMS_DIR INSERT_SIZE INTERVALBEDFILE REFFASTA CLAMMS_SPECIAL_REGIONS LIBRARY_DIR

#A faire tourner avec "chr", enlever "chr" en sortie du windows.bed
#le github préconise de faire tourner sans le "chr" => prevoir un genome fasta sans chr
#pour l'insert size prevoir un peu plus long que le fragment size ex: 200pb -> 250pb

windowsBed() {

  export CLAMMS_DIR=$1
  export INSERT_SIZE=$2
  INTERVALBEDFILE=$3
  REFFASTA=$4
  #REF = CLAMMS_DIR/lib4CLAMMS/hg19
  CLAMMS_SPECIAL_REGIONS=$5
  #CLAMMS_DIR/data
  LIBRARY_DIR=$6
  
  debug "windowsBed : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "windowsBed : INSERT_SIZE is : \"${INSERT_SIZE}\""
  debug "windowsBed : INTERVALBEDFILE is : \"${INTERVALBEDFILE}\""
  debug "windowsBed : REFFASTA is : \"${REFFASTA}\""
  debug "windowsBed : CLAMMS_SPECIAL_REGIONS is : \"${CLAMMS_SPECIAL_REGIONS}\""
  debug "windowsBed : LIBRARY_DIR is : \"${LIBRARY_DIR}\""

  info "Checking windowsBed's arguments..."
  # - If exists
  ## - Check if /PATH/TO/CLAMS exists
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR}\" is not a correct path to clams directory"
    help 
  fi 
  ## - Check if /PATH/TO/annotate_windows.sh exists 
  if [ ! -f "${CLAMMS_DIR}annotate_windows.sh" ]
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
  ## - Check if LIBRARY_DIR exists 
  if [ ! -d ${LIBRARY_DIR} ]
  then 
    error "\"${LIBRARY_DIR}\" does not exist !"
    help 
  fi

  # - Extension checking
  
  if [[ ${INTERVALBEDFILE##*\.} != "bed" ]]
  then 
    error "\"${INTERVALBEDFILE}\" must be .bed file !"
  fi 

  if [[ ${REFFASTA##*\.} != "fa" ]]
  then
    error "\"${REFFASTA}\" must be .fa file !"
  fi 

  if [[ ${CLAMMS_SPECIAL_REGIONS##*\.} != "bed" ]]
  then 
    error "\"${CLAMMS_SPECIAL_REGIONS}\" must be .bed file !"
  fi 
  info "... Argument checking : done !"

  info "Launching annotate_windows.sh ..."
  chmod +x ${CLAMMS_DIR}annotate_windows.sh
  # - Mkdir if insertSize${INSERT_SIZE} does not exist

if [ ! -d ${LIBRARY_DIR}windowsBeds ]
  then
    echo "\"${LIBRARY_DIR}windowsBeds \" does not exist but was created !  "
    mkdir ${LIBRARY_DIR}windowsBeds
  fi
  
if [ ! -d ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE} ]
  then 
    echo "\"${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}\" does not exist but was created !  "
    mkdir ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}
  fi 
  # - Sort INTERVALBEDFILE and sed chr column of the resulting file
  sort -k1,1 -k2,2n ${INTERVALBEDFILE} | sed 's/^chr//g'> ${LIBRARY_DIR}interval_sort_nochr.bed

  # - Run annotate_windows
  ${CLAMMS_DIR}annotate_windows.sh ${LIBRARY_DIR}interval_sort_nochr.bed ${REFFASTA} ${CLAMMS_DIR}/lib4Clamms/hg19/mappability.bed ${INSERT_SIZE} ${CLAMMS_SPECIAL_REGIONS} > ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}/windows_nochr_${INSERT_SIZE}pb.bed
  info "... annotate_windows.sh done !"
}



###########################################################
# NORMALIZE COVERAGE DATA
###########################################################

# -- After Samtools 

normalizeFS() {

  # - Command ./clams.sh normalizeFS /PATH/TO/coverage /PATH/TO/clams_dir /PATH/TO/WINDOWS_BED

  export COVERAGE_PATH=$1
  export CLAMMS_DIR=$2
  export WINDOWS_BED=$3
  export LIBRARY_DIR=$4

  debug "normalizeFS : COVERAGE_PATH is : \"${COVERAGE_PATH}\""
  debug "normalizeFS : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "normalizeFS : WINDOWS_BED is : \"${WINDOWS_BED}\""
  debug "normalizeFS : LIBRARY_DIR : \"${LIBRARY_DIR}\""
 
  info "Checking normalize's arguments ..."

  # - Check if Coverage Path exists 
  if [ ! -d ${COVERAGE_PATH} ]
  then
    error "\"${COVERAGE_PATH}\" does not exist !"
    help 
  fi 
  
  # - Check if CLAMMS_DIR exists
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR}\" does not exist ! "
    help 
  fi 

  # - Check if WINDOWS_BED exists 
  if [ ! -f ${WINDOWS_BED} ]
  then
    error "\"${WINDOWS_BED}\" does not exist ! "
    help
  fi

  info "... Argument checking : done !"
  info "Lauching normalize and removing "chr" in first column ..."

  cd ${COVERAGE_PATH}
  for i in  *.bed; do sed 's/^chr//g' $i > $(basename "$i" | cut -d_ -f1 | cut -d. -f1).clamms.coverage.bed ; done
  ls *.clamms.coverage.bed | cut -d '.' -f 1 | while read SAMPLEID ; do  ${CLAMMS_DIR}normalize_coverage ${SAMPLEID}.clamms.coverage.bed $WINDOWS_BED > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed ;done

  info "... normalize Done !"
}

normalize(){

  # - Command ./clams.sh normalize /PATH/TO/coverage /PATH/TO/clams_dir SAMPLE 

  CLAMMS_DIR=$1
  SAMPLEID=$2
  CLAMMSCOVERAGEFILE=$3
  WINDOWS_BED=$4
  LIBRARY_DIR=$5

  debug "normalize : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "normalize : SAMPLEID name is : \"${SAMPLEID}\""
  debug "normalize : CLAMMSCOVERAGEFILE is : \"${CLAMMSCOVERAGEFILE}\""
  debug "normalize : WINDOWS_BED is : \"${WINDOWS_BED}\""
  debug "normalize : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
 
  info "Checking normalize's arguments ..."
  # - If exists
  ## - Check if Coverage Path exists 
  if [ ! -d ${COVERAGE_PATH} ]
  then
    error "\"${COVERAGE_PATH}\" does not exist !"
    help 
  fi
  ## - Check 
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR}\" does not exist ! "
    help 
  fi 
  
  # - Extension checking 
  if [[ ${CLAMMSCOVERAGEFILE##*\.} != "bed" ]]
  then 
    error "\"${CLAMMSCOVERAGEFILE}\" must be .bed file !"
  fi 
 
  if [[ ${WINDOWS_BED##*\.} != "bed" ]]
  then
    error "\"${WINDOWS_BED}\" must be .bed file !"
  info "... Argument checking : done !"
  fi

  info "Lauching normalize ..."
  ${CLAMMS_DIR}normalize_coverage ${CLAMMSCOVERAGEFILE} ${WINDOWS_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed

  info "... normalize Done !"


}


###########################################
# PICARDMERGED - TO DO 
###########################################

#create metrics matrix for KD tree, all in one: (be careful, column position of picard metrics could change, better to do a grep instead)
metricsMatrixFS(){

  # - Command ./clams.sh metricsMatrix /PATH/TO/hs_metrics SAMPLE

  LIBRARY_DIR=$1
  HS_FOLDER=$2
  MATCH_METRICS=$3

  debug "metricsMatrixFS : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "metricsMatrixFS : HS_FOLDER is : \"${HS_FOLDER}\""

  info "Checking metricsMatrixFS' arguments ..."

  # - Check if HS_FOLDER exists 
  if [ ! -d ${HS_FOLDER} ]
  then 
    error "\"${HS_FOLDER}\" does not exist !"
    help 
  fi 

  info "... Argument checking : done !"
  info "Launching metricsMatrixFS ..."
  
  ## TOUTES LES VARIABLES SONT A RECUP DANS MOBIDL
  for i in ${HS_FOLDER}*hs_metrics.txt
  do 
    SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
    echo "SAMPLE" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt
    echo ${SAMPLEID} >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt 
    #path to multiple metrics
    grep -v "#" $i | head -3 | tail -2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_metrics.txt 
  done
  
for i in ${HS_FOLDER}*insertsize_metrics.txt
  do
    SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
    grep -v "#" $i | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_metrics.txt  >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_insertsize_metrics.txt
  done

    # Sample ID insertsizer metrics .txt | sample id hs metrics .txt 
    grep -v "#" $i | head -3 | tail -2 | cut -f5 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_insertsize_metrics.txt 
    paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_hsmetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_insertsize_metrics.txt >  ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt 

  #head the file with parameter names
  echo "SAMPLE	PCT_PF_UQ_READS	ON_BAIT_VS_SELECTED	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_50X	AT_DROPOUT	GC_DROPOUT	MEAN_INSERT_SIZE"  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt
  #fill with the data (check if conversion of "," into "." is nedded)
  for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*kdTree_metrics.txt; do tail -1 $i | sed 's/,/./g' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt ; done

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
  
  LIBRARY_DIR=$1
  SAMPLEID=$2
  HSMETRICSTXT=$3
  INSERT_SIZE_METRICS_TXT=$4

  debug "metricsMatrix : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "metricsMatrix : SAMPLEID is : \"${SAMPLEID}\""
  debug "metricsMatrix : HSMETRICSTXT is : \"${HSMETRICSTXT}\""
  debug "metricsMatrix : INSERT_SIZE_METRICS_TXT is : \"${INSERT_SIZE_METRICS_TXT}\""

  info "Checking metricsMatrix's arguments ..."


  info "... Argument checking : done !"
  info "Launching metricsMatrix ..."
  
  ## TOUTES LES VARIABLES SONT A RECUP DANS MOBIDL
  echo "SAMPLE" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt
  echo ${SAMPLEID} >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt 
  #path to multiple metrics
  grep -v "#" ${HSMETRICSTXT} | head -3 | tail -2 | cut -f51,42,38,21,10,52 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_hsmetrics.txt 
  #path to insertsize metrics
  # Sample ID insertsizer metrics .txt | sample id hs metrics .txt 
  grep -v "#" ${INSERT_SIZE_METRICS_TXT} | head -3 | tail -2 | cut -f5 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_insertsize_metrics.txt 
  paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_hsmetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_insertsize_metrics.txt >  ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt 
  # $SAMPLEID_kdTree_metrics.txt a ecrire dans le répertoire kdTreeMetrics
  
  #head the file with parameter names
  #cat ${SAMPLEID}_kdTree_metrics.txt  > ALL_kdTreeMetrics.txt
  #fill with the data (check if conversion of "," into "." is nedded)
  ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt t -1 $i | sed -i '1i' ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt


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

  RSCRIPT=$1
  RSCRIPT_FILE=$2
  KNN=$3 
  ALL_TREE=$4
  KD_OUT=$5
  FROM_SCRATCH=$6
  
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

  if [ ! -f ${RSCRIPT_FILE} ]
  then 
    error "\"${RSCRIPT_FILE}\" does not exist ! Please select a correct file."
    help 
  fi 

  if [[ ${RSCRIPT_FILE##*\.} != "rscript" ]]
  then 
    error "\"${RSCRIPT_FILE}\" is not a correct file ! Please select .rscript file."
    help 
  fi 

  if [ ! -f ${ALL_TREE} ]
  then 
    error "\"${ALL_TREE}\" does not exist ! Please select a correct file."
    help 
  fi 

  if [[ (${FROM_SCRATCH} != "fromscratch") && (${FROM_SCRATCH} != "") ]]
  then 
    error "\"${FROM_SCRATCH}\" is not correct. Please enter \"fromscratch\" if it is the first time you run this script or nothing if it is not."
    help 
  fi 

  info "... Argument checking : done !"

  info "Launching KdTree RScript ..."

  
  ${RSCRIPT} ${RSCRIPT_FILE} ${KNN} ${ALL_TREE} ${KD_OUT} ${FROM_SCRATCH}

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
  normalizeFS $2 $3 $4 $5
fi

if [ $1 == "normalize" ]
then 
  normalize $2 $3 $4 $5 $6
fi 

if [ $1 == "metricsMatrixFS" ]
then 
  metricsMatrixFS $2 $3
fi 

if [ $1 == "metricsMatrix" ]
then 
  metricsMatrix $2 $3 $4 $5 
fi 

if [ $1 == "dirpreparation" ]
then 
  case $2 in 
    clamms)
      dirpreparation clamms $3
      ;;
    library)
      dirpreparation library $3 $4
      ;;
    size)
      dirpreparation size $3 $4
      ;;
  esac
fi 
