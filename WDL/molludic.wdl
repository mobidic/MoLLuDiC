import "modules/CNVdirPreparation.wdl" as runCNVDirPreparation
import "modules/cnvFromScratch.wdl" as runCNVFromScratch
import "modules/normalize.wdl" as runNormalize
import "modules/metricsMatrix.wdl" as runMetricsMatrix
import "modules/removeRelatives.wdl" as runRemoveRelatives
import "modules/makekdtree.wdl" as runMakekdtree
import "modules/cnvCalling.wdl" as runCNVCalling
import "modules/annotationCNV.wdl" as runAnnotationCNV

workflow molludic {

  # Variables

  # Call

}
