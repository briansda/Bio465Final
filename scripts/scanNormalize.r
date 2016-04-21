library(SCAN.UPC)

sampleId = commandArgs()[7]
datasetId = commandArgs()[8]
probeBrainArrayPackage = commandArgs()[9]

inputFilePath = paste("/fslhome/aguyer/fsl_groups/fslg_Bio465_BDAG/geneExpression/", datasetId, "/unprocessed/", sampleId ,".CEL", sep = "")
scanOutputFilePath = paste("/fslhome/aguyer/fsl_groups/fslg_Bio465_BDAG/geneExpression/", datasetId, "/processed/", sampleId , sep = "")
#upcOutputFilePath = paste("/fslhome/aguyer/fsl_groups/fslg_piccololab/compute/Biomarker_Benchmark/GeneExpression/", datasetId, "/Processed/upc/upc_", sampleId ,".txt", sep = "")

print("SAMPLE: ")
print(sampleId)
print("DATASET: ")
print(datasetId)
print(scanOutputFilePath)
#print(upcOutputFilePath)
if (file.exists(scanOutputFilePath))
{
  message(paste("A SCAN output file already exists at ", scanOutputFilePath, ".", sep=""))
  stop()
}

#if (file.exists(upcOutputFilePath))
#{
 # message(paste("A UPC output file already exists at ", upcOutputFilePath, ".", sep=""))
 # stop()
#}

#library(doParallel)
#registerDoParallel(cores=11)

scanOutput = exprs(SCAN(inputFilePath, probeSummaryPackage = probeBrainArrayPackage))
#upcOutput = exprs(UPC(inputFilePath, probeSummaryPackage = probeBrainArrayPackage))

scanRowsToRemove = grepl("AFFX", rownames(scanOutput))
#upcRowsToRemove = grepl("AFFX", rownames(upcOutput))

scanOutput = scanOutput[which(!scanRowsToRemove),,drop=FALSE]
#upcOutput = upcOutput[which(!upcRowsToRemove),,drop=FALSE]

colnames(scanOutput) = sampleId
#colnames(upcOutput) = sampleId

rownames(scanOutput) = gsub("_at", "", rownames(scanOutput))
#rownames(upcOutput) = gsub("_at", "", rownames(upcOutput))

if (!dir.exists(dirname(scanOutputFilePath)))
  dir.create(dirname(scanOutputFilePath))

#if (!dir.exists(dirname(upcOutputFilePath)))
#  dir.create(dirname(upcOutputFilePath))

write.table(scanOutput, scanOutputFilePath, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
#write.table(upcOutput, upcOutputFilePath, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
