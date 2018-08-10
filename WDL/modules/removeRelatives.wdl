task removeRelatives {

  String SrunLow
  File ClammsWorkflow
  String LibraryDir
  File AllKdTree
  File FamilyList

  command {
  
    ${SrunLow} ${ClammsWorkflow} removeRelatives ${LibraryDir} ${AllKdTree} ${FamilyList}
  
  }

  output {

    

  }

}
