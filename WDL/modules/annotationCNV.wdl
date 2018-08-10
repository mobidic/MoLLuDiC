task annotationCNV {

  String SrunLow
  File ClammsWorkflow
  String LibraryDir
  String SampleID
  String BedToolsPath
  File HGBed
  File HeaderAnnotation
  File CNVBed

  command {

    ${SrunLow} ${ClammsWorkflow} annotation ${LibraryDir} ${SampleID} ${BedToolsPath} ${HGBed} ${HeaderAnnotation} ${CNVBed} 

  }

  output {

    File cnvAnnotated4Achab = "${LibraryDir}projects/all/normCoverageNoChrBeds/${SampleID}.annotated.forachab.bed"
    File cnvAnnotated = "${LibraryDir}projects/all/normCoverageNoChrBeds/${SampleID}.annotated.final.bed"

  }

}
