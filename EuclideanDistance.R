suppressMessages(library(dplyr))
library(readr)
library(MASS)

setwd("/Users/Brian/Desktop/data/")
geneExpression<- read_tsv('geneExpression2.txt')

matrix_thing <-read_tsv('data_4.txt')
rownames(matrix_thing)=matrix_thing[,1]
as.list(matrix_thing[,1])
as.character(matrix_thing[,1])
row.names(matrix_thing)
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
  print(files[1])
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


testNumP <- function(data, cutoff)
{
  hold = data < cutoff
  for(i in seq(2, ncol(hold))){
    print(sum(hold[,i]))
  }
}

selectTop <- function(subtype, colNum)
{
  if(subtype[99,colNum]==subtype[100,colNum])
  {
    counter=0
    for(i in seq(100,nrow(subtype)))
    {
      if(subtype[i,colNum]!=subtype[99,colNum])
      {
        counter=i
        break
      }
    }
    return(subtype[seq(1,counter-1),c(1,colNum)])
  }
  return(subtype[seq(1,100),c(1,colNum)])
}


pairwiseUnion <- function(cancer1, cancer2)
{
  output<-NULL
  #multi<-1
  hold <- inner_join(cancer1,cancer2)
  total=nrow(cancer1)+nrow(cancer2)
  div=total-nrow(hold)
  #if (nrow(cancer1)>nrow(cancer2))
  #{
  #  multi=nrow(cancer2)/nrow(cancer1)
  #}
  #else
  #{
  #  multi=nrow(cancer1)/nrow(cancer2)
  #}
  output<-c((nrow(hold)/div),nrow(hold))
#   for(i in seq(1,ncol(cancer1),by=2))
#   {
#     subtype1 = cancer1[,c(i,i+1)]
#     outputList<-NULL
#     for(j in seq(1,ncol(cancer2),by=2))
#     {
#       
#       subtype2 = cancer2[,c(j,j+1)]
#       hold <- inner_join(subtype1, subtype2)
#       if(nrow(hold) != 0)
#       {
#         outputList <- c(outputList,nrow(hold)*multi)
#       }
#     }
#     output <- cbind(output,outputList)
#   }
#   return(output)
}

#Breast
breastCancerPath <- "/Users/Brian/Desktop/data/21653_Breast/"
breast <- cancerTypeFactorSection(breastCancerPath, geneExpression)
hist(c(breast$X21653_Basal.txt,breast$X21653_Erbb2.txt,breast$X21653_LuminalA.txt,breast$X21653_LuminalB.txt), xlab="Range of p-values", main="Distribution of p-values from Breast Cancer Wilcoxon Tests")
breastA <- arrange(breast, X21653_Basal.txt)
breastA<-selectTop(breastA,2)
#breastA <- breastA[seq(1,1000),c(1,2)]

breastB <- arrange(breast, X21653_Erbb2.txt)
breastB<-selectTop(breastB,3)
#breastB <- breastB[seq(1,1000),c(1,3)]

breastC <- arrange(breast, X21653_LuminalA.txt)
breastC<-selectTop(breastC,4)
#breastC <- breastC[seq(1,1000),c(1,4)]

breastD <- arrange(breast, X21653_LuminalB.txt)
breastD<-selectTop(breastD,5)
#breastD <- breastD[seq(1,1000),c(1,5)]

#breastSelected <- cbind(breastA,breastB,breastC,breastD)


#Testis
testisCancerPath <- "/Users/Brian/Desktop/data/3218_Testis"
testis <- cancerTypeFactorSection(testisCancerPath, geneExpression)
hist(c(testis$X3218_EC.txt ,testis$X3218_mixed.txt,testis$X3218_Sem.txt, testis$X3218_Ter.txt , testis$X3218_YS.txt), xlab="Range of p-values", main="Distribution of p-values from Testicular Cancer Wilcoxon Tests")
testisA <- arrange(testis, X3218_AnapSem.txt)
testisA<-selectTop(testisA,2)
#testisA <- testisA[seq(1,1000),c(1,2)]

testisB <- arrange(testis, X3218_Chorio.txt)
testisB<-selectTop(testisB,3)
#testisB <- testisB[seq(1,1000),c(1,3)]

#testisC <- arrange(testis, X3218_EC.Sem.txt)
#testisC <- testisC[seq(1,100),c(1,4)]

#testisD <- arrange(testis, X3218_EC.Ter.txt)
#testisD <- testisD[seq(1,100),c(1,5)]

#testisE <- arrange(testis, X3218_EC.YS.txt)
#testisE <- testisE[seq(1,100),c(1,6)]

testisC <- arrange(testis, X3218_EC.txt)
testisC<-selectTop(testisC,4)
#testisC <- testisC[seq(1,1000),c(1,4)]

testisD <- arrange(testis, X3218_mixed.txt)
testisD<-selectTop(testisD,5)
#testisD <- testisD[seq(1,1000),c(1,5)]

testisE <- arrange(testis, X3218_Sem.txt)
testisE<-selectTop(testisE,6)
#testisE <- testisE[seq(1,1000),c(1,6)]

testisF <- arrange(testis, X3218_Ter.txt)
testisF<-selectTop(testisF,7)
#testisF <- testisF[seq(1,1000),c(1,7)]

testisG <- arrange(testis, X3218_YS.txt)
testisG<-selectTop(testisG,8)
#testisG <- testisG[seq(1,1000),c(1,8)]

#testisSelected <- cbind(testisA,testisB,testisC,testisD,testisE,testisF,testisG)

#ovarian
ovarianCancerPath <- "/Users/Brian/Desktop/data/6008_Ovarian"
ovarian <- cancerTypeFactorSection(ovarianCancerPath, geneExpression)
testis <- cancerTypeFactorSection(testisCancerPath, geneExpression)
hist(c(ovarian$X6008_ClearCell.txt, ovarian$X6008_Endometriod.txt, ovarian$X6008_Mucinous.txt, ovarian$X6008_Serous.txt), xlab="Range of p-values", main="Distribution of p-values from Ovarian Cancer Wilcoxon Tests")


ovarianA <- arrange(ovarian, X6008_ClearCell.txt)
ovarianA<-selectTop(ovarianA,2)
#ovarianA <- ovarianA[seq(1,1000),c(1,2)]

ovarianB <- arrange(ovarian, X6008_Endometriod.txt)
ovarianB<-selectTop(ovarianB,3)
#ovarianB <- ovarianB[seq(1,1000),c(1,3)]

ovarianC <- arrange(ovarian, X6008_Mucinous.txt)
ovarianC<-selectTop(ovarianC,4)
#ovarianC <- ovarianC[seq(1,1000), c(1,4)]

ovarianD <- arrange(ovarian, X6008_Serous.txt)
ovarianD<-selectTop(ovarianD,5)
#ovarianD <- ovarianD[seq(1,1000),c(1,5)]

#ovarianSelected <- cbind(ovarianA,ovarianB,ovarianC, ovarianD)

#Kidney
kidneyCancerPath <- "/Users/Brian/Desktop/data/7023_Kidney"
kidney <- cancerTypeFactorSection(kidneyCancerPath, geneExpression)
ovarian <- cancerTypeFactorSection(ovarianCancerPath, geneExpression)
testis <- cancerTypeFactorSection(testisCancerPath, geneExpression)
hist(c(kidney$X7023_Pap1.txt, kidney$X7023_Pap12.txt, kidney$X7023_Pap2A.txt, kidney$X7023_Pap2B.txt), xlab="Range of p-values", main="Distribution of p-values from Renal Cancer Wilcoxon Tests")


kidneyA <- arrange(kidney, X7023_Pap1.txt)
kidneyA<-selectTop(kidneyA,2)
#kidneyA <- kidneyA[,c(1,2)]

kidneyB <- arrange(kidney, X7023_Pap12.txt)
kidneyB<-selectTop(kidneyB,3)
#kidneyB <- kidneyB[,c(1,3)]

kidneyC <- arrange(kidney, X7023_Pap2A.txt)
kidneyC<-selectTop(kidneyC,4)
#kidneyC <- kidneyC[,c(1,4)]

kidneyD <- arrange(kidney, X7023_Pap2B.txt)
kidneyD<-selectTop(kidneyD,5)
#kidneyD <- kidneyD[,c(1,5)]

#kidneySelected <- cbind(kidneyA,kidneyB,kidneyC,kidneyD)

#Print out matrixes of data of pairwise comparisions of selected genes for each subtype within each cancer
#Best if ran a single line at a time as to better see output and not all at the same time.
test<-pairwiseUnion(breastA,breastA)
test
test<-pairwiseUnion(breastA,breastB)
test
test<-pairwiseUnion(breastA,breastC)
test
test<-pairwiseUnion(breastA,breastD)
test
test<-pairwiseUnion(breastB,breastC)
test
test<-pairwiseUnion(breastB,breastD)
test
test<-pairwiseUnion(breastC,breastD)
test

test<-pairwiseUnion(breastA,kidneyA)
test
test<-pairwiseUnion(breastA,kidneyB)
test
test<-pairwiseUnion(breastA,kidneyC)
test
test<-pairwiseUnion(breastA,kidneyD)
test
test<-pairwiseUnion(breastB,kidneyA)
test
test<-pairwiseUnion(breastB,kidneyB)
test
test<-pairwiseUnion(breastB,kidneyC)
test
test<-pairwiseUnion(breastB,kidneyD)
test
test<-pairwiseUnion(breastC,kidneyA)
test
test<-pairwiseUnion(breastC,kidneyB)
test
test<-pairwiseUnion(breastC,kidneyC)
test
test<-pairwiseUnion(breastC,kidneyD)
test
test<-pairwiseUnion(breastD,kidneyA)
test
test<-pairwiseUnion(breastD,kidneyB)
test
test<-pairwiseUnion(breastD,kidneyC)
test
test<-pairwiseUnion(breastD,kidneyD)
test

test<-pairwiseUnion(breastA,ovarianA)
test
test<-pairwiseUnion(breastA,ovarianB)
test
test<-pairwiseUnion(breastA,ovarianC)
test
test<-pairwiseUnion(breastA,ovarianD)
test
test<-pairwiseUnion(breastB,ovarianA)
test
test<-pairwiseUnion(breastB,ovarianB)
test
test<-pairwiseUnion(breastB,ovarianC)
test
test<-pairwiseUnion(breastB,ovarianD)
test
test<-pairwiseUnion(breastC,ovarianA)
test
test<-pairwiseUnion(breastC,ovarianB)
test
test<-pairwiseUnion(breastC,ovarianC)
test
test<-pairwiseUnion(breastC,ovarianD)
test
test<-pairwiseUnion(breastD,ovarianA)
test
test<-pairwiseUnion(breastD,ovarianB)
test
test<-pairwiseUnion(breastD,ovarianC)
test
test<-pairwiseUnion(breastD,ovarianD)
test

test<-pairwiseUnion(breastA,testisA)
test
test<-pairwiseUnion(breastA,testisB)
test
test<-pairwiseUnion(breastA,testisC)
test
test<-pairwiseUnion(breastA,testisD)
test
test<-pairwiseUnion(breastA,testisE)
test
test<-pairwiseUnion(breastA,testisF)
test
test<-pairwiseUnion(breastA,testisG)
test
test<-pairwiseUnion(breastB,testisA)
test
test<-pairwiseUnion(breastB,testisB)
test
test<-pairwiseUnion(breastB,testisC)
test
test<-pairwiseUnion(breastB,testisD)
test
test<-pairwiseUnion(breastB,testisE)
test
test<-pairwiseUnion(breastB,testisF)
test
test<-pairwiseUnion(breastB,testisG)
test
test<-pairwiseUnion(breastC,testisA)
test
test<-pairwiseUnion(breastC,testisB)
test
test<-pairwiseUnion(breastC,testisC)
test
test<-pairwiseUnion(breastC,testisD)
test
test<-pairwiseUnion(breastC,testisE)
test
test<-pairwiseUnion(breastC,testisF)
test
test<-pairwiseUnion(breastC,testisG)
test
test<-pairwiseUnion(breastD,testisA)
test
test<-pairwiseUnion(breastD,testisB)
test
test<-pairwiseUnion(breastD,testisC)
test
test<-pairwiseUnion(breastD,testisD)
test
test<-pairwiseUnion(breastD,testisE)
test
test<-pairwiseUnion(breastD,testisF)
test
test<-pairwiseUnion(breastD,testisG)
test


test<-pairwiseUnion(kidneyA,kidneyA)
test
test<-pairwiseUnion(kidneyA,kidneyB)
test
test<-pairwiseUnion(kidneyA,kidneyC)
test
test<-pairwiseUnion(kidneyA,kidneyD)
test
test<-pairwiseUnion(kidneyB,kidneyC)
test
test<-pairwiseUnion(kidneyB,kidneyD)
test
test<-pairwiseUnion(kidneyC,kidneyD)
test

test<-pairwiseUnion(kidneyA,ovarianA)
test
test<-pairwiseUnion(kidneyA,ovarianB)
test
test<-pairwiseUnion(kidneyA,ovarianC)
test
test<-pairwiseUnion(kidneyA,ovarianD)
test
test<-pairwiseUnion(kidneyB,ovarianA)
test
test<-pairwiseUnion(kidneyB,ovarianB)
test
test<-pairwiseUnion(kidneyB,ovarianC)
test
test<-pairwiseUnion(kidneyB,ovarianD)
test
test<-pairwiseUnion(kidneyC,ovarianA)
test
test<-pairwiseUnion(kidneyC,ovarianB)
test
test<-pairwiseUnion(kidneyC,ovarianC)
test
test<-pairwiseUnion(kidneyC,ovarianD)
test
test<-pairwiseUnion(kidneyD,ovarianA)
test
test<-pairwiseUnion(kidneyD,ovarianB)
test
test<-pairwiseUnion(kidneyD,ovarianC)
test
test<-pairwiseUnion(kidneyD,ovarianD)
test


test<-pairwiseUnion(kidneyA,testisC)
test
test<-pairwiseUnion(kidneyA,testisD)
test
test<-pairwiseUnion(kidneyA,testisE)
test
test<-pairwiseUnion(kidneyA,testisF)
test
test<-pairwiseUnion(kidneyA,testisG)
test
test<-pairwiseUnion(kidneyB,testisC)
test
test<-pairwiseUnion(kidneyB,testisD)
test
test<-pairwiseUnion(kidneyB,testisE)
test
test<-pairwiseUnion(kidneyB,testisF)
test
test<-pairwiseUnion(kidneyB,testisG)
test
test<-pairwiseUnion(kidneyC,testisC)
test
test<-pairwiseUnion(kidneyC,testisD)
test
test<-pairwiseUnion(kidneyC,testisE)
test
test<-pairwiseUnion(kidneyC,testisF)
test
test<-pairwiseUnion(kidneyC,testisG)
test
test<-pairwiseUnion(kidneyD,testisC)
test
test<-pairwiseUnion(kidneyD,testisD)
test
test<-pairwiseUnion(kidneyD,testisE)
test
test<-pairwiseUnion(kidneyD,testisF)
test
test<-pairwiseUnion(kidneyD,testisG)
test


test<-pairwiseUnion(ovarianA,ovarianA)
test
test<-pairwiseUnion(ovarianA,ovarianB)
test
test<-pairwiseUnion(ovarianA,ovarianC)
test
test<-pairwiseUnion(ovarianA,ovarianD)
test
test<-pairwiseUnion(ovarianB,ovarianC)
test
test<-pairwiseUnion(ovarianB,ovarianD)
test
test<-pairwiseUnion(ovarianC,ovarianD)
test


test<-pairwiseUnion(ovarianA,testisC)
test
test<-pairwiseUnion(ovarianA,testisD)
test
test<-pairwiseUnion(ovarianA,testisE)
test
test<-pairwiseUnion(ovarianA,testisF)
test
test<-pairwiseUnion(ovarianA,testisG)
test
test<-pairwiseUnion(ovarianB,testisC)
test
test<-pairwiseUnion(ovarianB,testisD)
test
test<-pairwiseUnion(ovarianB,testisE)
test
test<-pairwiseUnion(ovarianB,testisF)
test
test<-pairwiseUnion(ovarianB,testisG)
test
test<-pairwiseUnion(ovarianC,testisC)
test
test<-pairwiseUnion(ovarianC,testisD)
test
test<-pairwiseUnion(ovarianC,testisE)
test
test<-pairwiseUnion(ovarianC,testisF)
test
test<-pairwiseUnion(ovarianC,testisG)
test
test<-pairwiseUnion(ovarianD,testisC)
test
test<-pairwiseUnion(ovarianD,testisD)
test
test<-pairwiseUnion(ovarianD,testisE)
test
test<-pairwiseUnion(ovarianD,testisF)
test
test<-pairwiseUnion(ovarianD,testisG)
test


test<-pairwiseUnion(testisC,testisC)
test
test<-pairwiseUnion(testisC,testisD)
test
test<-pairwiseUnion(testisC,testisE)
test
test<-pairwiseUnion(testisC,testisF)
test
test<-pairwiseUnion(testisC,testisG)
test
test<-pairwiseUnion(testisD,testisD)
test
test<-pairwiseUnion(testisD,testisE)
test
test<-pairwiseUnion(testisD,testisF)
test
test<-pairwiseUnion(testisD,testisG)
test
test<-pairwiseUnion(testisE,testisE)
test
test<-pairwiseUnion(testisE,testisF)
test
test<-pairwiseUnion(testisE,testisG)
test
test<-pairwiseUnion(testisF,testisF)
test
test<-pairwiseUnion(testisF,testisG)
test
test<-pairwiseUnion(testisG,testisG)
test
