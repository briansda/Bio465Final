suppressMessages(library(dplyr))
library(readr)
library(MASS)

setwd("/Users/annaguyer/Desktop/data/")
geneExpression<- read_tsv('geneExpression.txt')

geneFactorSelection <- function(geneExpression, control, cancerSubtype)
{
  controlPatients <- read_tsv(control, col_names = FALSE)
  controlPatients <- as.vector(controlPatients[,1])
  cancerSubtypePatients <- read_tsv(cancerSubtype, col_names = FALSE)
  cancerSubtypePatients <- as.vector(cancerSubtypePatients[,1])
  pvals <- NULL
  for(i in seq(1, nrow(geneExpression)))
  {
    controlExp= NULL
    cancerExp= NULL
    for(controlP in controlPatients)
    {
      controlExp= c(controlExp, geneExpression[i,controlP])
    }
    for(cancerP in cancerSubtypePatients)
    {
      cancerExp= c(cancerExp, geneExpression[i,cancerP])
    }
    controlExp= as.numeric(controlExp)
    cancerExp= as.numeric(cancerExp)
    results <- wilcox.test(controlExp, cancerExp, exact= TRUE)
    pvals <- c(pvals, results$p.value)
  }
  return(pvals)
}

cancerTypeFactorSection <- function(path, geneExpression)
{
  pvalsMatrix <- geneExpression[,1]
  setwd(path)
  files <- list.files()
  fileNames = files[seq(2, length(files))]
  for(subFile in seq(2,length(files)))
  {
    control <- files[1]
    subtype <- files[subFile]
    #print(subtype)
    currPvals <- geneFactorSelection(geneExpression, control, subtype)
    #print(currPvals)
    if(is.null(pvalsMatrix))
    {
      pvalsMatrix = currPvals
    }else
    {
      pvalsMatrix = cbind(pvalsMatrix, currPvals)
    }
  }
  colnames(pvalsMatrix) = c("genes", fileNames)
  return(pvalsMatrix)
}

euclideanDistance <- function(geneExpression, patientsSubA, patientsSubB, genes)
{
  currGenes <- filter(geneExpression, Gene==genes)
  distr <- NULL
  for(patientSubA in patientsSubA)
  {
    patientSubAExp <- currGenes[,patientSubA]
    patientSubAExp <- as.vector(patientSubAExp)
    for(patientSubB in patientsSubB)
    {
      patientSubBExp <- currGenes[,patientSubB]
      patientSubBExp <- as.vector(patientSubBExp)
      hold <- dist(rbind(patientSubAExp, patientSubBExp), method = "euclidean")
      distr <- c(distr, hold)
    }
  }
  return(distr)
}


breastCancerPath <- "/Users/annaguyer/Desktop/data/21653_Breast/"
breast <- cancerTypeFactorSection(breastCancerPath, geneExpression)

testisCancerPath <- "/Users/annaguyer/Desktop/data/3218_Testis"
testis <- cancerTypeFactorSection(testisCancerPath, geneExpression)

ovarianCancerPath <- "/Users/annaguyer/Desktop/data/63885_Ovarian"
ovarian <- cancerTypeFactorSection(ovarianCancerPath, geneExpression)

kidneyCancerPath <- "/Users/annaguyer/Desktop/data/7023_Kidney"
kidney <- cancerTypeFactorSection(kidneyCancerPath, geneExpression)

testNumP <- function(data, cutoff)
{
  hold = data < cutoff
  for(i in seq(2, length(hold))){
    print(sum(hold[,i]))
  }
}
