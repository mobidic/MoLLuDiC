task dirpreparation {
  
  File ClammsWorkflow
  String Option
  String ClammsDir
  String LibraryName
  Float InsertSize 


  command {

    ${ClammsWorkflow} dirpreparation ${Option}${ClammsDir} ${LibraryName} ${InsertSize}

  }

}
