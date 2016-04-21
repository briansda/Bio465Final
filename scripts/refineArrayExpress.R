args <- commandArgs(TRUE)

refineData <- function(geo_id){
  filename = paste(geo_id,".txt",sep="")
  filename = paste("/Users/annaguyer/Desktop/UndergraduateResearch/AlgorithmAnalysis/ArrayExpress/", filename, sep="")
  geoDataset = read.delim(filename, header = TRUE, row.names=1, sep = '\t', quote="")
  geoDataset = geoDataset[grep("Characteristics|Factor.Value|Comment", names(geoDataset))]
  #View(geoDataset)
  toRemove = array()
  for (i in seq(1,length(names(geoDataset)))){
    #print(names(geoDataset)[i])
    names(geoDataset)[i] = stringr::str_replace(names(geoDataset)[i], pattern = "[:print:]+[.]{2}", replacement = "")
    names(geoDataset)[i] = stringr::str_replace(names(geoDataset)[i],"[.]", "")
    #check = as.character(unique(geoDataset[names(geoDataset)[i]]))
    #if (check=="1"){
    #  toRemove[length(toRemove)+1] = names(geoDataset)[i]
    #}
    check = grep(".file",names(geoDataset)[i])
    if (length(check) != 0){
      toRemove[length(toRemove)+1] = names(geoDataset)[i]
    }
  }
  geoDataset = geoDataset[, !duplicated(names(geoDataset))]
  geoDataset = geoDataset[,! names(geoDataset) %in% toRemove, drop = FALSE]
  filename = paste(geo_id,".tsv",sep="")
  filename = paste("/Users/annaguyer/Desktop/UndergraduateResearch/AlgorithmAnalysis/CleanArrayExpressClinicalData/", filename, sep="")
  
  write.table(geoDataset, file = filename, sep = "\t", col.names = NA, quote = FALSE, row.names = TRUE)

}
#print(args[1])
geo_ids <- strsplit(args[1], ",")[[1]]
#print(geo_ids)
for (i in (1:length(geo_ids))){
  geo_id = geo_ids[i]
  print(geo_id)
  print(i)
  geo_id = paste("GSE", geo_id, sep="")
  refineData(geo_id)
}

