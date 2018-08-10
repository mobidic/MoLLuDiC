task metricsMatrix {

  String SrunLow
  File ClammsWorkflow 
  String LibraryDir
  String SampleID
  File HSMetricsTxt
  File InsertSizeMetricsTxt #Bind this to MobiDL when available
  String PythonPath
  String KeepMetricsFiles

  command {

    ${SrunLow} ${ClammsWorkflow} metricsMatrix ${LibraryDir} ${SampleID} ${HSMetricsTxt} ${InsertSizeMetricsTxt} ${PythonPath} ${KeepMetricsFiles}

  }

  output {

    File allKdTree = "${LibraryDir}projects/all/kdTreeMetrics/ALL_kdTreeMetrics.txt"
    File kdTree ="${LibraryDir}projects/all/kdTreeMetrics/${SampleID}_ktTree_metrics.txt"

  }

}
