task normalizeFS {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  String LibraryName
  String CoveragePath
  File WindowsBedFile

  command {

    ${SrunLow} ${ClammsWorkflow} normalizeFS ${ClammsDir} ${LibraryName} ${CoveragePath} ${WindowsBedFile}

  }


}
