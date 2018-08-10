task windowsBed {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  String LibraryName
  Int InsertSize
  File IntervalBedFile
  File RefFasta
  File ClammsSpecialRegions

  command {
  
    ${SrunLow} ${ClammsWorkflow} windowsBed ${ClammsDir} ${LibraryName} ${InsertSize} ${IntervalBedFile} ${RefFasta} ${ClammsSpecialRegions}
  
  }

  output {
    
    File windowsBedFile = "${ClammsDir}/lib4Clamms/${LibraryName}/windowsBeds/insertSize${InsertSize}/windows_nochr_${InsertSize}pb.bed"

  }


}
