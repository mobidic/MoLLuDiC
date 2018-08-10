task makekdtree {

  String SrunLow
  String ClammsWorkflow
  String LibraryDir
  String SampleID
  String RscriptPath
  Int KNN
  File AllKdTreeRR

  command {
  
    ${SrunLow} ${ClammsWorkflow} makekdtree ${LibraryDir} ${RscriptPath} ${AllKdTreeRR}
  
  }

  output {

    File nnsFile = "${LibraryDir}projects/all/kdTreeMetrics/${SampleID}.${KNN}nns.txt"
      
  }

}
