#!/usr/bin/bash

###########################################################
#
# Bash script for clamms workflow 
#
# By MoBiDic - Version 0.0.2
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


VERSION="0.0.3"
VERBOSITY=4


###########################################################
# Help 
###########################################################

help() {
  echo "MoLLuDiC (version ${VERSION}) is a CNV workflow for calling and annotation !"
  echo "Usage : /.molludic.sh"
  echo "General arguments : "
  echo "    help : show this help message"
  echo "    -v : decrease of increase verbosity level (ERROR : 1 | WARNING : 2 | INFO [default] : 3 | DEBUG : 4)"
  echo "    "
  echo "MoLLuDiC is composed of several functions. You print help for each module by typing help after function name."
  echo "    Example : ./molludic.sh install help"a
  echo "List of MoLLuDiC's functions : "
  echo "  dirpreparation <OPTION> : Create folders to use correctly clamms"
  echo "  install <CLAMM_DIRECTORY> : install Clamms in specific directory"
  echo "  mapinstall <CLAMM_DIRECTORY> <BigWigToWig_PATH> : install Mapability bed"
  echo "  windowsBed <CLAMMS_DIRECTORY> <INSERT_SIZE> <INTERVALBEDFILE> <REFFASTA> <CLAMMS_SPECIAL_REGIONS> <LIBRARY_DIRECTORY> : run clamms annotate windows"
  echo "  normalizeFS <CLAMMS_DIRECTORY> <LIBRARY_NAME> <CLAMMSCOVERAGEFILE> <WINDOWS_BED> : normalize bed files from scratch"
  echo "  normalize <CLAMMS_DIRECTORY> <SAMPLEID> <CLAMMSCOVERAGEFILE> <WINDOWS_BED> <LIBRARY_DIRECTORY> : normalize one bed file"
  echo "  metricsMatrixFS <LIBRARY_DIRECTORY> <HS_FOLDER> <PYTHON_PATH> <MATCH_METRICS> : create kd tre metrics from scratch"
  echo "  metricsMatrix : <LIBRARY_DIRECTORY> <SAMPLEID> <HSMETRICSTXT> <INSERT_SIZE_METRICS_TXT> <PYTHON_PATH> <MATCH_METRICS> : create kd tree metric for 1 sample"
  echo "  removeRelatives <ALLKDTREE> <FAMILYLIST> <LIBRARY_DIRECTORY> : remove relatives from all kd tree file"
  echo "  makekdtree <RSCRIPT_PATH> <RSCRIPT_FILE> <KNN> <ALL_TREE> <LIBRARY_DIRECTORY> <FROM_SCRATCH> : use Rscript to do kd tree"
  echo "  cnvCallingFS <CLAMMS_DIRECTORY> <LIBRARY_DIRECTORY> <LIST_KDTREE> <WINDOWS_BED> <KNN> : do calling from scratch"
  echo "  cnvCalling <CLAMMS_DIRECTORY> <LIBRARY_DIRECTORY> <NORMCOVBED> <LIST_KDTREE> <WINDOWS_BED> <KNN> : do calling for 1 sample"
  echo "  annotation <LIBRARY_DIRECTORY> <SAMPLEID> <BEDTOOLS_PATH> <HGBED> <HEADER_FILE> <CNV_BED> <DAD> (optional) <MUM> (optiona) : annotate cnv bed file"
  echo " "
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
      if [ ! -d ${LIBRARY_DIR}windowsBeds ]
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

    all)
      CLAMMS_DIR=$2
      LIBRARY_NAME=$3
      INSER_SIZE=$4
      LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}
      ddebug "dirpreparation library : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
      debug "dirpreparation library : LIBRARY_NAME is : \"${LIBRARY_NAME}\""
      mkdir ${CLAMMS_DIR}
      mkdir ${CLAMMS_DIR}lib4Clamms
      mkdir ${CLAMMS_DIR}lib4Clamms/hg19
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
      if  [ ! -d ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds ]
      then 
        error "Invalid Library directory format : \"${LIBRARY_DIR}projects/all/normCoverageNoChrBeds\" is missing ! Please use : \"dirpreparation library CLAMMS_DIR LIBRARY_NAME\" to create correct library directory."
      fi
      if [ ! -d ${LIBRARY_DIR}windowsBeds ]
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

    help)
      echo "MoLLuDiC (version ${VERSION}) dirpreparation help !"
      echo "Usage : ./molludic.sh dirpreparation <OPTION> <other argument depending on OPTION."
      echo "Several OPTION available :"
      echo "    1 - ./molludic.sh dirpreparation clamms <CLAMMS_DIRECTORY> : With clamms option you will create library directories for clamms in <CLAMMS_DIRECTORY>. Use it if it is your first time."
      echo "    2 - ./molludic.sh dirpreparation library <CLAMMS_DIRECTORY> <LIBRARY_NAME> : With library option you will create your <LIBRARY_NAME> in your <CLAMMS_DIRECTORY>. Use it if you want to create a new library."
      echo "    3 - ./molludic.sh dirpreparation size <LIBRARY_DIRECTORY> <INSER_SIZE> : With size option you will create insert size folder in your <LIBRARY_DIRECTORY>, generated with library option (2). Use it if you want tu use new insert size."
      echo "    4 - ./molludic.sh dirpreparation all <CLAMMS_DIRECTORY> <LIBRARY_NAME> <INSERT_SIZE> : Do all previous steps. Use it if you want to create your CLAMMS_DIRECTORY, LIBRARY and INSERT_SIZE folder at once."
      echo "    5 - ./molludic.sh dirpreparation help : To print this message."
      echo " "
      exit 1
      ;;

    *)
      error "Invalid argument \"$1\" !"
      help
      ;;

  esac
}


###########################################################
# INSTALL 
###########################################################

install() {

  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) install help !"
    echo "Usage : ./molludic.sh install <CLAMMS_DIRECTORY>"
    echo "    <CLAMMS_DIRECTORY> : Directory where you want to install Clamms."
    echo " "
    exit 1
  else

    CLAMMS_DIR=$1
    debug "install : Installation PATH of clamms is : \"${CLAMMS_DIR}\""
    info "install : Checking installation directory of clamms..."

    # - Check if /PATH/TO/Install exists else create directory
    if [ ! -d ${CLAMMS_DIR} ] 
    then 
      warning "\"${CLAMMS_DIR}\" does not exist. \"${CLAMMS_DIR}\" was created"
      dirpreparation clamms ${CLAMMS_DIR}
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
  fi
}

###########################################################
# MAPABILITY INSTALL 
###########################################################

mapinstall() {

  if [ $1 =="help" ]
  then
    echo "MoLLuDiC (version ${VERSION} mapinstall help !"
    echo "Usage ./molludic.sh mapinstall <CLAMMS_DIRECTORY> <BW2W_PATH>"
    echo "    <CLAMMS_DIRECTORY> : Directory of Clamms."
    echo "    <BW2W_PATH> : Path to BigWigToWig."
    echo " "
    exit 1

  else 
  
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
      dirpreparation clamms ${CLAMMS_DIR}
    fi 
    # - Check if /PATH/TO/bigWigToWig exists
    if [ ! -f ${BW2W_PATH} ]
    then 
      error "\"${BW2W_PATH}\" is not a correct PATH to BigWigToWig !"
      mapinstall help 
    fi 
    # - Check if /PATH/TO/bigWigToWig is a PATH to the soft
    if ! [[ "${BW2W_PATH}" =~ bigWigToWig$ ]]
    then 
      error "\"${BW2W_PATH}\" is not bigWigToWig ! Please chose a correct PATH."
      mapinstall help 
    fi 
    info "... Arguments checking done"

    info "Installing Mapability ..."
    cd ${INSTALLATION_PATH}
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
    ${BW2W_PATH} wgEncodeCrgMapabilityAlign100mer.bigWig wgEncodeCrgMapabilityAlign100mer.wig
    grep -v '^#' wgEncodeCrgMapabilityAlign100mer.wig | sed 's/^chr//g' > mappability.bed
    info "... Done !"

  fi 

}


###########################################################
# CREATE WINDOWS BED
###########################################################

windowsBed() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) windowsBed help !"
    echo "Usage ./molludic.sh windowsbed <CLAMMS_DIRECTORY> <LIBRARY_NAME> <INSERT_SIZE> <INTERVALBEDFILE> <REFFASTA> <CLAMMS_SPECIAL_REGIONS>"
    echo "    <CLAMMS_DIRECTORY> : Directory where Clamms is installed."
    echo "    <LIBRARY_NAME> : CNV Library you want to use."
    echo "    <INSERT_SIZE> : Size of insert."
    echo "    <INTERVALBEDFILE> : Bed file of the library."
    echo "    <REFFASTA> : Reference genome in fasta format."
    echo "    <CLAMMS_SPECIAL_REGIONS> : Bed file in git repo of Clamms."
    exit 1
  else 

    export CLAMMS_DIR=$1
    LIBRARY_NAME=$2
    export INSERT_SIZE=$3
    INTERVALBEDFILE=$4
    REFFASTA=$5
    CLAMMS_SPECIAL_REGIONS=$6
    LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/

    debug "windowsBed : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
    debug "windowsBed : INSERT_SIZE is : \"${INSERT_SIZE}\""
    debug "windowsBed : INTERVALBEDFILE is : \"${INTERVALBEDFILE}\""
    debug "windowsBed : REFFASTA is : \"${REFFASTA}\""
    debug "windowsBed : CLAMMS_SPECIAL_REGIONS is : \"${CLAMMS_SPECIAL_REGIONS}\""
    debug "windowsBed : LIBRARY_NAME is : \"${LIBRARY_NAME}\""

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
    # - Mkdir if insertSize${INSERT_SIZE} does not exist

    if [ ! -d ${LIBRARY_DIR}windowsBeds ]
    then
      warning "\"${LIBRARY_DIR}windowsBeds \" does not exist but was created !  "
      mkdir ${LIBRARY_DIR}windowsBeds
    fi
    
    if [ ! -d ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE} ]
    then 
      warning "\"${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}\" does not exist but was created !  "
      mkdir ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}
    fi 
    # - Sort INTERVALBEDFILE and sed chr column of the resulting file
    sort -k1,1 -k2,2n ${INTERVALBEDFILE} | sed 's/^chr//g'> ${LIBRARY_DIR}interval_sort_nochr.bed

    # - Run annotate_windows
    ${CLAMMS_DIR}annotate_windows.sh ${LIBRARY_DIR}interval_sort_nochr.bed ${REFFASTA} ${CLAMMS_DIR}/lib4Clamms/hg19/mappability.bed ${INSERT_SIZE} ${CLAMMS_SPECIAL_REGIONS} > ${LIBRARY_DIR}windowsBeds/insertSize${INSERT_SIZE}/windows_nochr_${INSERT_SIZE}pb.bed
    info "... annotate_windows.sh done !"
  fi 
}



###########################################################
# NORMALIZE COVERAGE DATA
###########################################################


normalizeFS() {

  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) normalizeFS help !"
    echo "Usage ./molludic.sh normalizeFS <CLAMMS_DIRECTORY> <LIBRARY_NAME> <COVERAGE_PATH> <WINDOWS_BED>"
    echo "    <CLAMMS_DIRECTORY> : Directory where Clamms is installed."
    echo "    <LIBRARY_NAME> : CNV library you want to use."
    echo "    <COVERAGE_PATH> : Path to directory which contains coverage bed files."
    echo "    <WINDOWS_BED> : Bed file generated by windowsBed function (./molludic.sh windowsBed help for more informations)."
    exit 1
  else

    export CLAMMS_DIR=$1
    LIBRARY_NAME=$2
    export COVERAGE_PATH=$3
    export WINDOWS_BED=$4
    export LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}/

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
    # Think to do, remove or not theese files in COVERAGE_PATH
    for i in ${COVERAGE_PATH}*.bed; do sed 's/^chr//g' $i > ${COVERAGE_PATH}$(basename "$i" | cut -d_ -f1 | cut -d. -f1).clamms.coverage.bed ; done
    for i in ${COVERAGE_PATH}*.clamms.coverage.bed ; do  ${CLAMMS_DIR}normalize_coverage ${COVERAGE_PATH}$(basename "$i" | cut -d_ -f1 | cut -d. -f1).clamms.coverage.bed $WINDOWS_BED > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/$(basename "$i" | cut -d_ -f1 | cut -d. -f1).norm.coverage.bed ;done

    info "... normalize Done !"
  fi 
}

normalize(){

  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) normalizeFS help !"
    echo "Usage ./molludic.sh normalize <CLAMMS_DIRECTORY> <LIBRARY_NAME> <SAMPLEID> <CLAMMSCOVERAGEFILE> <WINDOWS_BED>"
    echo "    <CLAMMS_DIRECTORY> : Directory where Clamms is installed."
    echo "    <LIBRARY_NAME> : CNV library you want to use."
    echo "    <SAMPLEID> : ID of wanted case"
    echo "    <CLAMMSCOVERAGEFILE> : Coverage bed file of patient."
    echo "    <WINDOWS_BED> : Bed file generated by windowsBed function (./molludic.sh windowsBed help for more informations)."
    exit 1
  else


    CLAMMS_DIR=$1
    LIBRARY_NAME=$2
    SAMPLEID=$3
    CLAMMSCOVERAGEFILE=$4
    WINDOWS_BED=$5
    LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}

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
  fi 


}

###########################################
# PICARDMERGED
###########################################

metricsMatrixFS(){

  if [ $1 == "help" ]
  then
    echo "MoLLuDiC (version ${VERSION}) metricsMatrixFS help !"
    echo "Usage ./molludic.sh metricsMatrixFS <LIBRARY_DIRECTORY> <HS_FOLDER> <PYTHON_PATH> <KEEP_FILE>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <HS_FOLDER> : Path to your folder containing HS metric files."
    echo "    <PYTHON_PATH> : Path to your current Python version (Python3 recommanded). You can use python3 if binded."
    echo "    <KEEP_FILE> (OPTIONNAL) : KEEP if you want to keep all files, \"\" else."
    exit 1
  else 

    LIBRARY_DIR=$1
    HS_FOLDER=$2
    PYTHON_PATH=$3

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
      grep "Y" ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "100"; else print "1"; }' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt
      
      # Create a compatible hsmetrics file
      grep -v "#" $i | head -3 | tail -2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt 
      # add insertsizemetrics
      grep -v "#" ${HS_FOLDER}${SAMPLEID}*insertsize_metrics.txt | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt -  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt
      info "Metrics are matching for ${SAMPLEID}"
     # select only useful column
      ${PYTHON_PATH} matchkdmetrics.py ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt | head -n2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt
     # add samplei and sex to metrics
     paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sample.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kdTree_metrics.txt 
    done
    #head the file with parameter names
    echo "SAMPLE	PCT_PF_UQ_READS	ON_BAIT_VS_SELECTED	PCT_TARGET_BASES_10X	PCT_TARGET_BASES_50X	AT_DROPOUT	GC_DROPOUT	MEAN_INSERT_SIZE	SEX"  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt
    #fill with the data (check if conversion of "," into "." is nedded)
    for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*kdTree_metrics.txt; do tail -n +2 $i | sed 's/,/./g' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt ; done
    info "... metricsMatrixFS done !"

    if [ "$4" != "KEEP" ]
    then
     info "remove temporary files"
      rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_tmp_*
      rm ${LIBRARY_DIR}projects/all/kdTreeMetrics/*_kd_sample.txt
    fi
  fi

}

metricsMatrix() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) metricsMatrix help !"
    echo "Usage ./molludic.sh metricsMatrix <LIBRARY_DIRECTORY> <SAMPLEID> <HSMETRICSTXT> <INSERT_SIZE_METRICS_TXT> <PYTHON_PATH> <KEEP>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <SAMPLEID> : ID of wanted case."
    echo "    <HSMETRICSTXT> : Metrics file for your case."
    echo "    <INSERT_SIZE_METRICS_TXT> : File generated by MobiDL containing your insert size metrics."
    echo "    <PYTHON_PATH> : Path to your current Python version (Python3 recommanded). You can use python3 if binded."
    echo "    <KEEP_FILE> (OPTIONNAL) : KEEP if you want to keep all files, \"\" else."
    exit 1
  else

    LIBRARY_DIR=$1
    SAMPLEID=$2
    HSMETRICSTXT=$3
    INSERT_SIZE_METRICS_TXT=$4
    PYTHON_PATH=$5

    debug "metricsMatrix : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
    debug "metricsMatrix : SAMPLEID is : \"${SAMPLEID}\""
    debug "metricsMatrix : HSMETRICSTXT is : \"${HSMETRICSTXT}\""
    debug "metricsMatrix : INSERT_SIZE_METRICS_TXT is : \"${INSERT_SIZE_METRICS_TXT}\""


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
    # - Create kd_sex_sample.txt
    echo "SEX" > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt
    grep "Y" ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.norm.coverage.bed | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "100"; else print "1"; }' >> ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_kd_sex_sample.txt

    #path to multiple metrics
    grep -v "#" ${HSMETRICSTXT} | head -3 | tail -2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt

    #path to insertsize metrics
    grep -v "#" ${INSERT_SIZE_METRICS_TXT} | head -3 | tail -2  | paste ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_hsmetrics.txt -  > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt

    info "Metrics are matching for ${SAMPLEID}"

    # select only useful column
     ${PYTHON_PATH} matchkdmetrics.py ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_metrics.txt | head -n2 > ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLEID}_tmp_kd_allmetrics.txt
 
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
  fi 

}


##########################################
# REMOVE RELATIVES 
##########################################

removeRelatives() {

  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) removeRelatives help !"
    echo "Usage ./molludic.sh removeRelatives <LIBRARY_DIRECTORY> <ALLKDTREE> <FAMILYLIST>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <ALLKDTREE> : ALL_kdTreeMetrics.txt generated by metricsMatrix(FS) functions (./molludic.sh metricsMatrixFS help for more informations)."
    echo "    <FAMILYLIST> : Members of case's family."
    exit 1
  else 

    export  LIBRARY_DIR=$1
    export  ALLKDTREE=$2
    export  FAMILYLIST=$3

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

    info "removeRelatives done !"
  fi 
}

###########################################
# MAKEKDTREE
###########################################

makekdtreeFS() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) makekdtree help !"
    echo "Usage ./molludic.sh makekdtreeFS <LIBRARY_DIRECTORY> <RSCRIPT_PATH> <KNN>" #<ALLKDTREE> <FROM_SCRATCH>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <RSCRIPT_PATH> : Path to your current Rscript version. You can use Rscript if binded."
    echo "    <KNN> : Your KNN's value."
    #echo "    <ALLKDTREE> : ALL_kdTreeMetrics.txt generated by metricsMatrix(FS) functions (./molludic.sh metricsMatrixFS help for more informations)."
    #echo "    <FROM_SCRATCH> (OPTIONNAL) : \"fromscratch\" if you want to keep all files, \"\" else."
    exit 1
  else

    LIBRARY_DIR=$1
    RSCRIPT_PATH=$2
    KNN=$3 
    #ALLKDTREE=$4
    #FROM_SCRATCH=$5
    
    debug "makekdtreeFS : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
    debug "makekdtreeFS : RSCRIPT_PATH is : \"${RSCRIPT_PATH}\""
    debug "makekdtreeFS : KNN is : \"${KNN}\""
    #debug "makekdtree : ALLKDTREE is : \"${ALLKDTREE}\""
    #debug "makekdtree : FROM_SCRATCH is : \"${FROM_SCRATCH}\""

    info "Checking makekdtree's arguments ..."

    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      makekdtreeFS help 
    fi 

    if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
    then 
      error "\"${KNN}\" is not a correct value ! Must be an integer."
      makekdtreeFS help  
    fi 

    #if [ ! -f ${ALLKDTREE} ]
    #then 
    #  error "\"${ALLKDTREE}\" does not exist ! Please select a correct file."
    #  makekdtree help
    #fi 

    #if [[ (${FROM_SCRATCH} != "fromscratch") && (${FROM_SCRATCH} != "") ]]
    #then 
    #  error "\"${FROM_SCRATCH}\" is not correct. Please enter \"fromscratch\" if it is the first time you run this script or nothing if it is not."
    #  makekdtree help
    #fi 

    info "... Argument checking : done !"

    info "Launching KdTreeFS RScript ..."

 
    for SAMPLE in ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLE}_ALL_kdTreeMetrics.txt
    do
      ${RSCRIPT_PATH} compute_kdTree.R ${KNN} ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLE}_ALL_kdTreeMetrics.txt ${LIBRARY_DIR}projects/all/kdTreeMetrics/ #${FROM_SCRATCH}
    done 

    #sort nns.txt for the next JOIN
    for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*.${KNN}nns.txt
    do 
      sort $i > ${i/txt/sort.txt}
      rm $i
    done

    info "... KdTreeFS RScript done !"
  fi 
}


makekdtree() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) makekdtree help !"
    echo "Usage ./molludic.sh makekdtree <LIBRARY_DIRECTORY> <RSCRIPT_PATH> <KNN> <ALLKDTREE>" #<FROM_SCRATCH>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <RSCRIPT_PATH> : Path to your current Rscript version. You can use Rscript if binded."
    echo "    <KNN> : Your KNN's value."
    echo "    <ALLKDTREE> : ALL_kdTreeMetrics.txt generated by metricsMatrix(FS) functions (./molludic.sh metricsMatrixFS help for more informations)."
    #echo "    <FROM_SCRATCH> (OPTIONNAL) : \"fromscratch\" if you want to keep all files, \"\" else."
    exit 1
  else

    LIBRARY_DIR=$1
    RSCRIPT_PATH=$2
    KNN=$3 
    ALLKDTREE=$4
    #FROM_SCRATCH=$5
    
    debug "makekdtree : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
    debug "makekdtree : RSCRIPT_PATH is : \"${RSCRIPT_PATH}\""
    debug "makekdtree : KNN is : \"${KNN}\""
    debug "makekdtree : ALLKDTREE is : \"${ALLKDTREE}\""
    #debug "makekdtree : FROM_SCRATCH is : \"${FROM_SCRATCH}\""

    info "Checking makekdtree's arguments ..."

    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      makekdtree help 
    fi 

    if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
    then 
      error "\"${KNN}\" is not a correct value ! Must be an integer."
      makekdtree help  
    fi 

    if [ ! -f ${ALLKDTREE} ]
    then 
      error "\"${ALLKDTREE}\" does not exist ! Please select a correct file."
      makekdtree help
    fi 

    #if [[ (${FROM_SCRATCH} != "fromscratch") && (${FROM_SCRATCH} != "") ]]
    #then 
    #  error "\"${FROM_SCRATCH}\" is not correct. Please enter \"fromscratch\" if it is the first time you run this script or nothing if it is not."
    #  makekdtree help
    #fi 

    info "... Argument checking : done !"

    info "Launching KdTree RScript ..."

 
    ${RSCRIPT_PATH} compute_kdTree.R ${KNN} ${ALLKDTREE} ${LIBRARY_DIR}projects/all/kdTreeMetrics/ ${FROM_SCRATCH}

    #sort nns.txt for the next JOIN
    for i in ${LIBRARY_DIR}projects/all/kdTreeMetrics/*.${KNN}nns.txt
    do 
      sort $i > ${i/txt/sort.txt}
      rm $i
    done

    info "... KdTree RScript done !"
  fi 
}

############################################
# MODEL + CNV CALLING 
############################################

cnvCallingFS() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) cnvCallingFS help !"
    echo "Usage ./molludic.sh cnvCallingFS <CLAMMS_DIRECTORY> <LIBRARY_NAME> <LIST_KDTREE> <WINDOWS_BED> <KNN>"
    echo "    <CLAMMS_DIRECTORY> : Directory where Clamms is installed."
    echo "    <LIBRARY_NAME> : CNV Library you want to use."
    echo "    <WINDOWS_BED> : Bed file generated by windowsBed function (./molludic.sh windowsBed help for more informations)."
    echo "    <KNN> : Your KNN's value."
    exit 1
  else 

    CLAMMS_DIR=$1
    LIBRARY_NAME=$2
    LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}
    WINDOWS_BED=$3
    KNN=$4

    debug "cnvCalling : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
    debug "cnvCalling : LIBRARY_NALE is : \"${LIBRARY_NAME}"
    debug "cnvCalling : LIBRARY_DIR is :\"${LIBRARY_DIR}\""
    debug "cnvCalling : WINDOWS_BED is : \"${WINDOWS_BED}\""
    debug "cnvCalling : KNN is : \"${KNN}\""

    info "cnvCalling : Arguments checking ..."

    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      cnvCallingFS help
    fi 

    if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
    then 
      error "\"${KNN}\" is not a correct value ! Must be an integer."
      cnvCallingFS help
    fi 
    
    if [ ! -d ${CLAMMS_DIR} ]
    then 
      error "\"${CLAMMS_DIR} does not exist !"
      cnvCallingFS help
    fi 

    if [ ! -f ${WINDOWS_BED} ]
    then 
      error "\"${WINDOWS_BED}\" does not exist !"
      cnvCallingFS help
    fi 

    info "cnvCalling : Arguments are OK !"
    info "cnvCalling : Begin calling ..."


    ls ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/*.norm.coverage.bed |sort | while read FILE
    do 
      SAMPLE=$(basename "${FILE}" | cut -d_ -f1 | cut -d. -f1)
      echo -e -n "${SAMPLE}\t${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.norm.coverage.bed\t"
      grep "Y" ${FILE} | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'
    done > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt

    ls ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/*.norm.coverage.bed | cut -d '.' -f 1 | while read FILE
    do 
      SAMPLE=$(basename "${FILE}" | cut -d_ -f1 | cut -d. -f1)
      SEX=`echo "${SAMPLE}" | join - ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 3`
      join ${LIBRARY_DIR}projects/all/kdTreeMetrics/${SAMPLE}.${KNN}nns.sort.txt ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 2- > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt
      ${CLAMMS_DIR}fit_models ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt ${WINDOWS_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed
      ${CLAMMS_DIR}call_cnv ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.norm.coverage.bed ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed --sex ${SEX} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.cnv.bed 
    done
    

    info "cnvCalling : Calling done !"
  fi 

}


cnvCalling() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) cnvCallingFS help !"
    echo "Usage ./molludic.sh cnvCalling <CLAMMS_DIRECTORY> <LIBRARY_NAME> <NORMCOVBED> <WINDOWS_BED> <NNS_PATH> <KNN>"
    echo "    <CLAMMS_DIRECTORY> : Directory where Clamms is installed."
    echo "    <LIBRARY_NAME> : CNV Library you want to use."
    echo "    <NORMCOVBED> : Bed file generated by normalize(FS) functions (./molludic.sh normalizeFS help for more information)"
    echo "    <WINDOWS_BED> : Bed file generated by windowsBed function (./molludic.sh windowsBed help for more informations)."
    echo "    <NNS_PATH> : nns.txt file generated by makekdtree function (./molludic.sh makekdtree help for more informations)."
    echo "    <KNN> : Your KNN's value."
    exit 1

    CLAMMS_DIR=$1
    LIBRARY_NAME=$2
    LIBRARY_DIR=${CLAMMS_DIR}lib4Clamms/${LIBRARY_NAME}
    NORMCOVBED=$3
    NNS_PATH=$4
    WINDOWS_BED=$5
    KNN=$6

    debug "cnvCalling : CLAMMS_DIR is : \"${CLAMMS_DIR}\""
    debug "cnvCalling : LIBRARY_DIR is :\"${LIBRARY_DIR}\""
    debug "cnvCalling : SAMPLEID is \"${SAMPLEID}\""
    debug "cnvCalling : NORMCOVBED is : \"${NORMCOVBED}\""
    debug "cnvCalling : NNS_PATH is : \"${NNS_PATH}\""
    debug "cnvCalling : WINDOWS_BED is : \"${WINDOWS_BED}\""
    debug "cnvCalling : KNN is : \"${KNN}\""

    info "cnvCalling : Arguments checking ..."

    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      cnvCalling help
    fi 

    if ! [[ "${KNN}" =~ ^[0-9]+$ ]]
    then 
      error "\"${KNN}\" is not a correct value ! Must be an integer." 
      cnvCalling help
    fi 
    
    if [ ! -d ${CLAMMS_DIR} ]
    then 
      error "\"${CLAMMS_DIR} does not exist !"
      cnvCalling help
    fi 

    if [ ! -f ${NNS_PATH} ]
    then 
      error "\"${NNS_PATH}\"  does not exist !"
      cnvCalling help
    fi 

    if [ ! -f ${WINDOWS_BED} ]
    then 
      error "\"${WINDOWS_BED}\" does not exist !"
      cnvCalling help
    fi 

    if [ ! -f ${NORMCOVBED} ]
    then 
      error "\"${NORMCOVBED}\" does not exist !"
      cnvCalling help
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
    join ${NNS_PATH} ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/all.file.sex.sort.txt | tr ' ' '\t' | cut -f 2- > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt
    ${CLAMMS_DIR}fit_models ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.${KNN}nns.ref.panel.sex.txt ${WINDOWS_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed
    ${CLAMMS_DIR}call_cnv ${NORMCOVBED} ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.models.bed --sex ${SEX} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLE}.cnv.bed 
  fi 

}


############################################
# ANNOTATION 
############################################

annotationFS() {

  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) annotationFS help !"
    echo "Usage ./molludic.sh annotationFS <LIBRARY_DIRECTORY> <BEDTOOLS_PATH> <HGBED> <CNV_BED>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <BEDTOOLS_PATH> : Path to Bedtools on your computer. You can use bedtools if binded."
    echo "    <HGBED> : HG19 annotated bed file."
    echo "    <CNV_BED> : cnv.bed file generated thanks to cnvCalling(FS) function (./molludic.sh cnvCallingFS help for more informations)."
    exit 1
  else 

    LIBRARY_DIR=$1
    BEDTOOLS_PATH=$2
    HGBED=$3 #/home/puce/resources/hg19/gene_annot_hg19_final.bed annotation.bed
    CNV_BED=$4
    
    info "annotationFS : argument checking ..."
    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      annotationFS help
    fi 

    if [ ! -f ${BEDTOOLS_PATH} ]
    then 
      error "\"${BEDTOOLS_PATH}\" does not exist !"
      annotationFS help 
    fi 

    if [ ! -f ${HGBED} ]
    then 
      error "\"${HGBED}\" does not exist !"
      annotationFS help 
    fi 

    if [ ! -f ${CNV_BED} ]
    then 
      error "\"${CNV_BED}\" does not exist !"
      annotationFS help
    fi 
    info "annotationFS : argument checking done !"
    debug "annotation : LIBRARY_DIR is : \"${LIBRARY_DIR}\""
    debug "annotation : BEDTOOLS_PATH is : \"${BEDTOOLS_PATH}\""
    debug "annotation : HGBED is : \"${HGBED}\""
    debug "annotation : CNV_BED is : \"${CNV_BED}\""

    info "Starting annotation ..."

    for i in ${CNV_BED}*.cnv.bed
    do 
      SAMPLEID=$(basename "$i" | cut -d_ -f1 | cut -d. -f1)
      info "Annotate for ${SAMPLEID}"
      # remove unwanted column 
      cut -f1-10 $i > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cnv.cut.bed
      # create a header
      cat headers/header.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed 
      cat headers/header.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.final.bed
      cat headers/header.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.filtered.final.bed

      #annotate with all data
      ${BEDTOOLS_PATH} intersect -a ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cnv.cut.bed -b ${HGBED} -loj >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed

      #add size of CNV + link to decipher genome browser
      awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}'  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.final.bed
      awk '$1 <= 50000 && $10 <= 10 {print $0}' ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.final.bed >> ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.filtered.final.bed
      awk -F "\t" '{print $16"\t"$8}' ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.filtered.final.bed > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.forachab.bed

    done
    info "Annotation done !"
  fi 

  }


annotation() {
  
  if [ $1 == "help" ]
  then 
    echo "MoLLuDiC (version ${VERSION}) annotationFS help !"
    echo "Usage ./molludic.sh annotationFS <LIBRARY_DIRECTORY> <BEDTOOLS_PATH> <HGBED> <HEADER_FILE> <CNV_BED> <DAD> <MUM>"
    echo "    <LIBRARY_DIRECTORY> : Path to your CNV library."
    echo "    <SAMPLEID> : ID of wanted case."
    echo "    <BEDTOOLS_PATH> : Path to Bedtools on your computer. You can use bedtools if binded."
    echo "    <HGBED> : HG19 annotated bed file."
    echo "    <HEADER_FILE> : Bed file in \"headers\" folder in MoLLuDiC git repo. Take header.herit.annotated.bed for trio, else header.annotated.bed."
    echo "    <CNV_BED> : cnv.bed file generated thanks to cnvCalling(FS) function (./molludic.sh cnvCallingFS help for more informations)."
    echo "    <DAD> (OPTIONAL) : DAD ID."
    echo "    <MUM> (OPTIONAL) : MUM ID."
    exit 1
  else 

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
    
    info "annotation : argument checking ..."
    if [ ! -d ${LIBRARY_DIR} ]
    then 
      error "\"${LIBRARY_DIR}\" does not exist !"
      annotation help
    fi 

    if [ ! -f ${BEDTOOLS_PATH} ]
    then 
      error "\"${BEDTOOLS_PATH}\" does not exist !"
      annotation help 
    fi 

    if [ ! -f ${HGBED} ]
    then 
      error "\"${HGBED}\" does not exist !"
      annotation help 
    fi 

    if [ ! -f ${HEADER_FILE} ]
    then 
      error "\"${HEADER_FILE}\" does not exist ! You can find these files in \"headers\" folder in MoLLuDiC git repo."
      annotation help 
    fi 

    if [ ! -f ${CNV_BED} ]
    then 
      error "\"${CNV_BED}\" does not exist !"
      annotation help
    fi 
    info "annotation : argument checking done !"


    info "Starting annotation ..."

    if [[ (${DAD} == "") || (${MUM} == "") ]]
    then
      debug "Not in Dad and Mum condition"
      # remove unwanted column 
      cut -f1-10 ${CNV_BED} > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cnv.cut.bed

      # create a header
      cat ${HEADER_FILE} >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed
      #annotate with all data
      ${BEDTOOLS_PATH} intersect -a ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.cnv.cut.bed -b ${HGBED} -loj >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed
      awk -F "\t" '{print $16"\t"$8}' ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.forachab.bed

      #add size of CNV + link to decipher genome browser
      awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}'  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.annotated.final.bed 
    else 
      debug "In Dad and Mum condition"
      ${BEDTOOLS_PATH} intersect -a ${CNV_BED} -b ${FAMILY_BED} -loj | sort -k1,1 -k2,2n | ${BEDTOOLS_PATH} merge -c 4,5,6,7,8,9,10,24,23,25 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse,collapse -delim "|" >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.bed 
      cat ${HEADER_FILE} >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed 
      ${BEDTOOLS_PATH} intersect -a  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.bed -b ${HGBED} -loj >>  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed
      awk -F "\t" '{print $19"\t"$8}' ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed > ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.forachab.bed

      awk 'BEGIN { FS = OFS = "\t" }{if(NR>1){print $3-$2,"https://decipher.sanger.ac.uk/browser#q/"$4,$0}else{print}}'  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.bed >  ${LIBRARY_DIR}projects/all/normCoverageNoChrBeds/${SAMPLEID}.HERIT-DN.annotated.final.bed 
    fi 

    info "Annotation done !"
  fi 

  }


###########################################################
# MAIN 
###########################################################

case $1 in 
  help)
    help
    ;;
  dirpreparation)
    case $2 in 
      clamms)
        dirpreation clamms $3
        ;;
      library)
        dirpreparation library $3 $4
        ;;
      size)
        dirpreparation size $3 $4
        ;;
      all)
        dirpreparation all $3 $4 $5
        ;;
      help)
        dirpreparation help
        ;;
      *)
        error "Incorrect arguments \"$2\""
        dirpreparation help;
        ;;
    esac
    ;;
  install)
    install $2
    ;;
  mapinstall)
    mapinstall $2 $3
    ;;
  windowsBed)
    windowsBed $2 $3 $4 $5 $6 $7
    ;;
  normalizeFS)
    normalizeFS $2 $3 $4 $5
    ;;
  normalize)
    normalize $2 $3 $4 $5 $6
    ;;
  metricsMatrixFS)
    metricsMatrixFS $2 $3 $4 $5
    ;;
  metricsMatrix)
    metricsMatrix $2 $3 $4 $5 $6 $7
    ;;
  removeRelatives)
    removeRelatives $2 $3 $4
    ;;
  makekdtreeFS)
    makekdtreeFS $2 $3 $4
    ;;
  makekdtree)
    makekdtree $2 $3 $4 $5
    ;;
  cnvCallingFS)
    cnvCallingFS $2 $3 $4 $5
    ;;
  cnvCalling)
    cnvCalling $2 $3 $4 $5 $6 $7
    ;;
  annotationFS)
    annotationFS $2 $3 $4 $5
    ;;
  annotation)
    annotation $2 $3 $4 $5 $6 $7 $8 $9
    ;;
esac 
