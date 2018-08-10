task metricsMatrixFS {
    
    String SrunLow
    File ClammsWorkflow 
    String LibraryDir
    String HSFolder
    String PythonPath
    String KeepMetricsFiles

    command {

      ${SrunLow} ${ClammsWorkflow} metricsMatrixFS ${LibraryDir} ${HSFolder} ${KeepMetricsFiles}

    }


}
