task cnvFromScratch {

  String SrunLow
  File ClammsWorkflow
  String ClammsDir
  String LibraryName
  Int InserSize
  File BW2WPath
  File IntervalBedFile
  File RefFasta
  File ClammsSpecialRegions
  String CoveragePath
  File WindowsBedFile
  String HSFolder
  String PythonPath
  String KeepMetricsFiles
  File AllKdTree
  File FamilyList
  String RscriptPath
  Int KNN
  File AllKdTreeRR
  String BedToolsPath
  File HGBed


  command {

    ${SrunLow} ${ClammsWorkflow} dirpreparation all${ClammsDir} ${LibraryName} ${InsertSize}
    ${ClammsWorkflow} install ${ClammsDir}
    ${ClammsWorkflow} mapinstall ${ClammsDir} ${BW2WPath}
    ${ClammsWorkflow} windowsBed ${ClammsDir} ${LibraryName} ${InsertSize} ${IntervalBedFile} ${RefFasta} ${ClammsSpecialRegions}
    ${ClammsWorkflow} normalizeFS ${ClammsDir} ${LibraryName} ${CoveragePath} ${WindowsBedFile}
    ${ClammsWorkflow} metricsMatrixFS ${LibraryDir} ${HSFolder} ${KeepMetricsFiles}
    ${ClammsWorkflow} removeRelatives ${LibraryDir} ${AllKdTree} ${FamilyList}
    ${ClammsWorkflow} makekdtreeFS ${LibraryDir} ${RscriptPath} ${KNN}
    ${ClammsWorkflow} cnvCallingFS ${ClammsDir} ${LibraryName} ${WindowsBedFile} ${KNN}
    ${ClammsWorkflow} annotationCNVFS ${LibraryDir} ${BedToolsPath} ${HGBed} ${LibraryDir}projects/all/normCoverageNoChrBeds/
    
  }


}
