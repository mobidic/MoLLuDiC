task install {

  File ClammsWorkflow
  String ClammsDir 

  command {
    
    ${ClammsWorkflow} install ${ClammsDir}
  
  }

}
