task mapinstall {
  
  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  File BW2WPath

  command {

    ${SrunLow} ${ClammsWorkflow} mapinstall ${ClammsDir} ${BW2WPath}

  }

}
