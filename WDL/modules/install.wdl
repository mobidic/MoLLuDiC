task install {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir 

  command {
    
    ${SrunLow} ${ClammsWorkflow} install ${ClammsDir}
  
  }

}
