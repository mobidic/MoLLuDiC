task annotationCNVTrio {

  String SrunLow
  File ClammsWorkflow
  String LibraryDir
  String SampleID
  String BedToolsPath
  File HGBed
  File HeaderAnnotation
  File CNVBed
  String FatherSample
  String MotherSample

  command {

    ${SrunLow} ${ClammsWorkflow} annotation ${LibraryDir} ${SampleID} ${BedToolsPath} ${HGBed} ${HeaderAnnotation} ${CNVBed} ${FatherSample} ${MotherSample}

  }

  output {

    File cnvAnnotatedTrio4Achab = "${LibraryDir}projects/all/normCoverageNoChrBeds/${SampleID}.HERIT-DN.annotated.forachab.bed"
    File cnvAnnotatedTrio = "${LibraryDir}projects/all/normCoverageNoChrBeds/${SampleID}.HERIT-DN.annotated.final.bed"

  }

}
