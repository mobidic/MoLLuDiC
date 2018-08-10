task normalize {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  String LibraryName
  String SampleID
  File ClammsCoverageFile
  File WindowsBedFile

  command {

    ${SrunLow} ${ClammsWorkflow} normalize ${ClammsDir} ${LibraryName} ${SampleID} ${ClammsCoverageFile} ${WindowsBedFile}

  }

  output {

    File normCovBed = "${ClammsDir}/lib4Clamms/${LibraryName}/projects/all/normCoverageNoChrBeds/${SampleID}.norm.coverage.bed"

  }

}
