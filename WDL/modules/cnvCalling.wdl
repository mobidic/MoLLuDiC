task cnvCalling {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  String LibraryName
  String SampleID
  File NormCovBed
  File NNSPath
  File WindowsBed
  Int KNN

  command {

    ${SrunLow} ${ClammsWorkflow} cnvCalling ${ClammsDir} ${LibraryName} ${NormCovBed} ${NNSPath} ${WindowsBed} ${KNN}

  }

  output {

    File cnvBed = "${ClammsDir}lib4Clamms/${LibraryName}projects/all/normCoverageNoChrBeds/${SampleID}.cnv.bed"

  }


}
