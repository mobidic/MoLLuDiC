task dirpreparation {
 
  String StrunLow
  File ClammsWorkflow
  String DirOption
  String ClammsDir
  String LibraryName
  Float InsertSize 


  command {

    ${SrunLow} ${ClammsWorkflow} dirpreparation ${DirOption}${ClammsDir} ${LibraryName} ${InsertSize}

  }

}
