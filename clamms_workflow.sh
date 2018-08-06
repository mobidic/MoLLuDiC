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

  # - Command ./clams.sh metricsMatrix /PATH/TO/LIBRARY_DIR /PATH/TO/HS_FOLDER /PATH/TO/MATCH_METRICS_SCRIPT 
  LIBRARY_DIR=$1
  HS_FOLDER=$2
  PYTHON_PATH=$3
  MATCH_METRICS=$4

  debug "metricsMatrixFS : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "metricsMatrixFS : HS_FOLDER is : \"${HS_FOLDER}\""
  debug "metricsMatrixFS : MATCH_METRICS_SCRIPT is : \"${MATCH_METRICS}\""

  info "Checking metricsMatrixFS' arguments ..."

  # - Check if HS_FOLDER exists 
  if [ ! -d ${HS_FOLDER} ]
  then 
    error "\"${HS_FOLDER}\" does not exist !"
    help 
  fi 

  if [ ! -d ${LIBRARY_DIR} ]
  then
    error "\"${LIBRARY_DIR}\" does not exist !"
    help
  fi

  if [ ! -f ${MATCH_METRICS} ]
  then
    error "\"${MATCH_METRICS}\" does not exist !"
    help
  fi

  info "... Argument checking : done !"
  info "Launching metricsMatrixFS ..."


 
 # - Create a combine metrics file for kdTree
  for i in ${HS_FOLDER}*hs_metrics.txt
  do 
    SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
    # Create kd_sample.txt
    echo "SAMPLE" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt
    echo ${SAMPLEID} >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt 
    # Create kd_sex_sample.txt
    echo "SEX" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt
    grep "Y" ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "100"; else print "0"; }' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt
    
    # Create a compatible hsmetrics file
    grep -v "#" $i | head -3 | tail -2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt 
    # add insertsizemetrics
    grep -v "#" ${HS_FOLDER}${SAMPLEID}*insertsize_metrics.txt | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt -  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt
 info "Metrics are matching for ${SAMPLEID}"
   # select only useful column
    ${PYTHON_PATH} ${MATCH_METRICS} ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt | head -n2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt
   # add samplei and sex to metrics
   paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt 
  done
  
 # for i in ${HS_FOLDER}*insertsize_metrics.txt
 # do
    #SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
   # grep -v "#" $i | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt -  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt
  #done

    # Sample ID insertsizer metrics .txt | sample id hs metrics .txt 
#  for i in ${HS_FOLDER}*tmp_kd_metrics.txt
#  do
#    SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
 #   python3 ${MATCH_METRICS} ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt
 # done

#  for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*kd_sample.txt
#  do 
#   SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1) 
#   paste $i ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt 
#  done
  # $SAMPLEID_kdTree_metrics.txt a ecrire dans le répertoire kdTreeMetrics

  #head the file with parameter names
  echo "SAMPLE	PCT_PF_UQ_READS	ON_BAIT_VS_SELECTED	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_50X	AT_DROPOUT	GC_DROPOUT	MEAN_INSERT_SIZE	SEX"  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt
  #fill with the data (check if conversion of "," into "." is nedded)
  for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*kdTree_metrics.txt; do tail -n +2 $i | sed 's/,/./g' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt ; done

  info "... metricsMatrixFS done !"

  if [ "$4" != "DEBUG" ]
  then
   info "remove temporary files"
  rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_tmp_*
  rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_kd_sample.txt
  fi

}

metricsMatrix() {
  
  # - Command ./clamms_workflow.sh metricsMatrix /PATH/TO/LIBRARY_DIR SAMPLEID /PATH/TO/HSMETRICSTXT /PATH/TO/INSERT_SIZE_METRICS_TXT /PATH/TO/MATCH_METRICS_SCRIPT 
  
  LIBRARY_DIR=$1
  SAMPLEID=$2
  HSMETRICSTXT=$3
  INSERT_SIZE_METRICS_TXT=$4
  PYTHON_PATH=$5
  MATCH_METRICS=$6

  debug "metricsMatrix : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "metricsMatrix : SAMPLEID is : \"${SAMPLEID}\""
  debug "metricsMatrix : HSMETRICSTXT is : \"${HSMETRICSTXT}\""
  debug "metricsMatrix : INSERT_SIZE_METRICS_TXT is : \"${INSERT_SIZE_METRICS_TXT}\""
  debug "metricsMatrix : MATCH_METRICS_SCRIPT is : \"${MATCH_METRICS}\""


  info "Checking metricsMatrix's arguments ..."

  if [ ! -f ${HSMETRICSTXT} ]
  then
    error "\"${HSMETRICSTXT}\" does not exist !"
    help
  fi

  if [ ! -f ${INSERT_SIZE_METRICS_TXT} ]
  then
    error "\"${INSERT_SIZE_METRICS_TXT}\" does not exist !"
    help
  fi

  if [ ! -d ${LIBRARY_DIR} ]
  then
    error "\"${LIBRARY_DIR}\" does not exist !"
    help
  fi

  if [ ! -f ${MATCH_METRICS} ]
  then
    error "\"${MATCH_METRICS}\" does not exist !"
    help
  fi

  info "... Argument checking : done !"
  info "Launching metricsMatrix ..."
  
 # - Create a combine metrics file for kdTree
  echo "SAMPLE" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt
  echo ${SAMPLEID} >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt 
 # Create kd_sex_sample.txt
  echo "SEX" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt
  grep "Y" ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "100"; else print "0"; }' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt

  #path to multiple metrics
  grep -v "#" ${HSMETRICSTXT} | head -3 | tail -2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt

  #path to insertsize metrics
  grep -v "#" ${INSERT_SIZE_METRICS_TXT} | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt -  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt

  info "Metrics are matching for ${SAMPLEID}"

  # select only useful column
   ${PYTHON_PATH} ${MATCH_METRICS} ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt | head -n2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt
 
   # add sample to metrics
   paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt

   # check if sample line already exist, if yes remove it
  sed -i "/^${SAMPLEID}\t/d" ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt  

  #fill with the data (check if conversion of "," into "." is nedded)
  tail -n +2 ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt | sed 's/,/./g' | sed -i '1 r /dev/stdin' ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt 

  info "... metricsMatrix done !"

  if [ "$6" != "DEBUG" ]
  then
   info "remove temporary files"
  rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_tmp_*
  rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_kd_sample.txt
  fi

}


##########################################
# REMOVE RELATIVES 
##########################################

removeRelatives() {

export  ALLKDTREE=$1
export  FAMILYLIST=$2
export  LIBRARY_DIR=$3

  debug "removeRelatives : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "removeRelatives : FAMILYLIST is : \"${FAMILYLIST}\""
  debug "removeRelatives : ALLKDTREE is : \"${ALLKDTREE}\""


  info "Checking removeRelatives's arguments ..."

  if [ ! -f ${ALLKDTREE} ]
  then
    error "\"${ALLKDTREE}\" does not exist !"
    help
  fi

  if [ ! -f ${FAMILYLIST} ]
  then
    error "\"${FAMILYLIST}\" does not exist !"
    help
  fi

  if [ ! -s ${FAMILYLIST} ]
  then
    error "\"${FAMILYLIST}\" is empty !"
  help
  exit 1
  fi

  if [ ! -d ${LIBRARY_DIR} ]
  then
    error "\"${LIBRARY_DIR}\" does not exist !"
    help
  fi

  info "... Argument checking : done !"
  info "Launching removeRelatives ..."

        awk 'BEGIN{family="";}{
                split($0,a,"\t");
                if(NF == 0 ){
                        next;
                }
                if(length(a)== 1){
                        system("cp ${ALLKDTREE} ${LIBRARY_DIR}projects/all/kdTreeMetrics/"a[1]"_ALL_kdTreeMetrics.txt" )
                }else{
                        for (i=1;i<=length(a);i++){
                                for (j=1;j<=length(a);j++){
                                        if(j!=i){family  = family"^"a[j]"\\s|"}
                                }
                                family = "\""family"^####################\"";
                                print a[i] "family is here:  "family;
                                system("grep -vE "family" ${ALLKDTREE} >  ${LIBRARY_DIR}projects/all/kdTreeMetrics/"a[i]"_ALL_kdTreeMetrics.txt")
                                family=""
                        };
                }

        }'  ${FAMILYLIST}
}

###########################################
# MAKEKDTREE
###########################################

makekdtree(){

  RSCRIPT_PATH=$1
  RSCRIPT_FILE=$2
  KNN=$3 
  ALL_TREE=$4
  LIBRARY_DIR=$5
  FROM_SCRATCH=$6
  
  debug "makekdtree : RSCRIPT_PATH is : \"${RSCRIPT_PATH}\""
  debug "makekdtree : RSCRIPT_FILE is : \"${RSCRIPT_FILE}\"" 
  debug "makekdtree : KNN is : \"${KNN}\""
  debug "makekdtree : ALL_TREE is : \"${ALL_TREE}\""
  debug "makekdtree : KD_OUT is : \"${LIBRARY_DIR}\""
  debug "makekdtree : FROM_SCRATCH is : \"${FROM_SCRATCH}\""

  info "Checking makekdtree's arguments ..."

  if [ ! -d ${LIBRARY_DIR} ]
  then 
    error "\"${LIBRARY_DIR}\" does not exist !"
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

  if [[ ${RSCRIPT_FILE##*\.} != "Rscript" ]]
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

 
 
  ${RSCRIPT_PATH} ${RSCRIPT_FILE} ${KNN} ${ALL_TREE} ${LIBRARY_DIR}projects/all/kdTreeMetrics/ ${FROM_SCRATCH}

  #sort nns.txt for the next JOIN
  for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*.${KNN}nns.txt
  do 
    sort $i > ${i/txt/sort.txt}
    rm $i
  done

  info "... KdTree RScript done !"
}


############################################
# MODEL + CNV CALLING 
############################################

cnvCallingFS() {
  
  CLAMMS_DIR=${1}
  LIBRARY_DIR=${2}
  LIST_KDTREE=${3}
  WINDOWS_BED=${4}
  KNN=${5}

  debug "cnvCalling : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "cnvCalling : LIBRARY_DIR is :\"${LIBRARY_DIR}\""
  debug "cnvCalling : LIST_KDTREE is : \"${LIST_KDTREE}\""
  debug "cnvCalling : WINDOWS_BED is : \"${WINDOWS_BED}\""
  debug "cnvCalling : KNN is : \"${KNN}\""

  info "cnvCalling : Arguments checking ..."

  if [ ! -d ${LIBRARY_DIR} ]
  then 
    error "\"${LIBRARY_DIR}\" does not exist !"
    help
  fi 

  if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
  then 
    error "\"${KNN}\" is not a correct value ! Must be an integer."
    help 
  fi 
  
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR} does not exist !"
    help 
  fi 

  if [ ! -f ${LIST_KDTREE} ]
  then 
    error "\"${LIST_KDTREE}\"  does not exist !"
    help 
  fi 

  if [ ! -f ${WINDOWS_BED} ]
  then 
    error "\"${WINDOWS_BED}\" does not exist !"
    help
  fi 

  info "cnvCalling : Arguments are OK !"
  info "cnvCalling : Begin calling ..."


  ls ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/*.norm.coverage.bed |sort | while read FILE
  do 
    SAMPLE=$(basename "${FILE}" | cut -d_ -f1 | cut -d. -f1)
    echo -e -n "${SAMPLE}\t${SAMPLE}.norm.coverage.bed\t"
    grep "Y" ${FILE} | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'
  done > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt

  ls ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/*.norm.coverage.bed | cut -d '.' -f 1 | while read FILE
  do 
    SAMPLE=$(basename "${FILE}" | cut -d_ -f1 | cut -d. -f1)
    SEX=`echo "${SAMPLE}" | join - ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 3`
    join ${LIST_KDTREE} ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 2- > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt
    ${CLAMMS_DIR}fit_models ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt ${WINDOWS_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed
    ${CLAMMS_DIR}call_cnv ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.norm.coverage.bed ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed --sex ${SEX} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.cnv.bed 
  done
  

  info "cnvCalling : Calling done !"

}


cnvCalling() {

  CLAMMS_DIR=${1}
  LIBRARY_DIR=${2}
  NORMCOVBED=${3}
  LIST_KDTREE=${4}
  WINDOWS_BED=${5}
  KNN=${6}

  debug "cnvCalling : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
  debug "cnvCalling : LIBRARY_DIR is :\"${LIBRARY_DIR}\""
  debug "cnvCalling : SAMPLEID is \"${SAMPLEID}\""
  debug "cnvCalling : NORMCOVBED is : \"${NORMCOVBED}\""
  debug "cnvCalling : LIST_KDTREE is : \"${LIST_KDTREE}\""
  debug "cnvCalling : WINDOWS_BED is : \"${WINDOWS_BED}\""
  debug "cnvCalling : KNN is : \"${KNN}\""

  info "cnvCalling : Arguments checking ..."

  if [ ! -d ${LIBRARY_DIR} ]
  then 
    error "\"${LIBRARY_DIR}\" does not exist !"
    help
  fi 

  if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
  then 
    error "\"${KNN}\" is not a correct value ! Must be an integer."
    help 
  fi 
  
  if [ ! -d ${CLAMMS_DIR} ]
  then 
    error "\"${CLAMMS_DIR} does not exist !"
    help 
  fi 

  if [ ! -f ${LIST_KDTREE} ]
  then 
    error "\"${LIST_KDTREE}\"  does not exist !"
    help 
  fi 

  if [ ! -f ${WINDOWS_BED} ]
  then 
    error "\"${WINDOWS_BED}\" does not exist !"
    help
  fi 

  if [ ! -f ${NORMCOVBED} ]
  then 
    error "\"${NORMCOVBED}\" does not exist !"
  fi 

  info "cnvCalling : Arguments are OK !"
  info "cnvCalling : Begin calling ..." 
  
  ls ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/*.norm.coverage.bed |sort | while read FILE
  do 
    SAMPLE=$(basename "${FILE}" | cut -d_ -f1 | cut -d. -f1)
    echo -e -n "${SAMPLE}\t${SAMPLE}.norm.coverage.bed\t"
    #${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${FILE}\t"
    grep "Y" ${FILE} | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'
  done > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt


  SAMPLE=$(basename "${NORMCOVBED}" | cut -d_ -f1 | cut -d. -f1)
  SEX=`echo "${SAMPLE}" | join - ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 3`
  join ${LIST_KDTREE} ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 2- > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt
  ${CLAMMS_DIR}fit_models ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt ${WINDOWS_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed
  ${CLAMMS_DIR}call_cnv ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.norm.coverage.bed ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed --sex ${SEX} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.cnv.bed 

}


############################################
# ANNOTATION 
############################################

annotation(){
  
  LIBRARY_DIR=$1
  export SAMPLEID=$2 #Sample_name_du_CI
  BEDTOOLS_PATH=$3
  HGBED=$4 #/home/puce/resources/hg19/gene_annot_hg19_final.bed annotation.bed
  HEADER_FILE=$5 #~/PROJECTS/EXOMES/header.annotated.bed or header.herit.annotated.bed 
  CNV_BED=$6
  export DAD=$7 #Sample_name_du_pere
  export MUM=$8 #Sample_name_du_mere

  debug "annotation : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
  debug "annotation : SAMPLEID is : \"${SAMPLEID}\""
  debug "annotation : BEDTOOLS_PATH is : \"${BEDTOOLS_PATH}\""
  debug "annotation : HGBED is : \"${HGBED}\""
  debug "annotation : HEADER_FILE is : \"${HEADER_FILE}\""
  debug "annotation : CNV_BED is : \"${CNV_BED}\""
  debug "annotation : DAD is : \"${DAD}\""
  debug "annotation : MUM is : \"${MUM}\""

  export FAMILY_LIST="${DAD} ${MUM}"   #Argument-liste à donner en debut de run 
  export FAMILY_BED=`sed 's/ /.cnv.bed /g; s/$/.cnv.bed/' <<< ${FAMILY_LIST}` #prevoir le chemin vers les bed
export FAMILY_GREP=`sed 's/ /|/g' <<< ${FAMILY_LIST}`

  debug "annotation : FAMILY_LIST is : \"${FAMILY_LIST}}\""
  debug "annotation : FAMILY_BED is : \"${FAMILY_BED}\""
  debug "annotation : FAMILY_GREP is \"${FAMILY_GREP}\""

  info "Starting annotation ..."

  if [[ (${DAD} == "") || (${MUM} == "") ]]
  then
    debug "Not in Dad and Mum condition"
    # remove unwanted column 
    cut -f1-10 ${CNV_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cut.cnv.bed

    # create a header
    cat ${HEADER_FILE} >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed
    #annotate with all data
    ${BEDTOOLS_PATH} intersect -a ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cut.cnv.bed -b ${HGBED} -loj >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed

    #add size of CNV + link to decipher genome browser
    awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}'  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.final.bed 
  else 
    debug "In Dad and Mum condition"
    ${BEDTOOLS_PATH} intersect -a ${CNV_BED} -b ${FAMILY_BED} -loj | sort -k1,1 -k2,2n | ${BEDTOOLS_PATH} merge -c 4,5,6,7,8,9,10,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse -delim "|" >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.bed 
    cat ${HEADER_FILE} >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed 
    ${BEDTOOLS_PATH} intersect -a  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.bed -b ${HGBED} -loj >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed

    awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}'  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.final.bed 
  fi 

  info "Annotation done !"

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
  metricsMatrixFS $2 $3 $4 $5
fi 

if [ $1 == "metricsMatrix" ]
then 
  metricsMatrix $2 $3 $4 $5 $6 $7
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

if [ $1 == "removeRelatives" ]
then
 removeRelatives $2 $3 $4
fi

if [ $1 == "makekdtree" ]
then 
  makekdtree $2 $3 $4 $5 $6 $7
fi

if [ ${1} == "cnvCallingFS" ]
then 
  cnvCallingFS ${2} ${3} ${4} ${5} ${6}
fi 

if [ $1 == "cnvCalling" ]
then 
  cnvCalling $2 $3 $4 $5 $6 $7 
fi 

if [ $1 == "annotation" ]
then 
  annotation $2 $3 $4 $5 $6 $7 $8 $9
fi 
