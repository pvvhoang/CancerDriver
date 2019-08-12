#================================================================
#' This function allows you to convert miRNA names to version 21.
#' @param miRs miRNA names
#' @return Converted miRNA names
#================================================================
convertmiRs <- function(miRs) {
  # Get version of miRs
  version = checkMiRNAVersion(miRs, verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature = miRNA_PrecursorToMature(miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miRNA list
  # If a converted item is "NA", get the corresponding value from the mature list
  miRs <- miRNameList_Converted[, 2]
  naList <- which(is.na(miRNameList_Converted[, 2]))
  for(i in 1:length(naList)){
    k <- naList[i]
    miRs[k] <- miRMature[k, 2]
  }
  
  return(miRs)
}

#================================================================
#' miRNA target prediction with the Pearson correlation coefficient method, returning p-value
#' Calculate the Pearson correlation coefficient of each pair of miRNA-mRNA,and return a matrix of p-value of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' 
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @return A  matrix that includes the p-value of Pearson correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' 
#' @references
#' Pearson, K. (1920) Notes on the history of correlation. Biometrika, 13, 25 - 45.
#================================================================
Pearson_pValue=function(datacsv, cause, effect){
  data=Read(datacsv)
  data=scale(data) #standardise the data
  
  header<-readHeader(datacsv)
  num_miRNA<-length(cause)
  num_mRNA<-length(effect)
  miR<-header[1:num_miRNA]
  mR<-header[-(1:num_miRNA)]
  
  miRNA=data[,cause]
  mRNA=data[,effect]
  
  r <- matrix(data = NA, nrow = num_mRNA, ncol = num_miRNA)
  row.names(r) <- mR
  colnames(r) <- miR
  for (i in 1:num_mRNA) {
    for (j in 1:num_miRNA) {
      r[i,j] <- cor.test(mRNA[,i], miRNA[,j], method="pearson")$p.value
    }
  }
  
  return(r)
}

#================================================================
#' Adjust p-values
#================================================================
adjustpValues = function(results) {
  
  r <- results
  
  nR <- nrow(r)
  nC <- ncol(r)
  t <- as.vector(r)
  t <- p.adjust(t, method="fdr")
  r <- matrix(t, nrow = nR, ncol = nC)
  
  row.names(r) <- row.names(results)
  colnames(r) <- colnames(results)
  
  return(r)
}

#================================================================
#' Identify links among nodes from their expression data
#' @param data Expression data with rows being samples and columns being biological features
#' @param cause Range of cause
#' @param effect Range of effect
#' @param rootDir Root folder
#' @param f File name
#' @param usingpVal TRUE if using p-value
#' @param cutoff FDR cutoff
#' @return Edges from cause to effect
#================================================================
identifyEdges = function(data, cause, effect, rootDir, f = "temp", usingpVal = FALSE, cutoff = 0.05) {
  
  if(usingpVal == TRUE) {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson_pValue(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Get links which have p-value < cutoff
    results <- adjustpValues(results)
    results <- results < cutoff
    ind <- which(results, arr.ind=TRUE)
    edges <- ind
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  } else {
    # Use Pearson to evaluate the relationship among cause and effect
    dataset <- paste(rootDir, "/Data/Output/dataset.csv", sep = "")
    write.csv(data[, c(cause, effect)], dataset, row.names = FALSE)
    results = Pearson(dataset, 1:length(cause), (length(cause)+1):(length(cause)+length(effect)))
    results <- t(results)
    
    # Write file
    t <- paste(rootDir, "/Data/Output/CancerDriver/Cancer/Network/", f, ".csv", sep = "")
    write.csv(results, t, row.names = TRUE)
    
    # Get links which have the absolute coefficients more than
    # the average of all absolute coefficients
    a <- mean(abs(results))
    results <- abs(results) >= a
    ind <- which(results, arr.ind=TRUE)
    edges <- ind
    edges[, 1] <- row.names(ind)
    edges[, 2] <- colnames(results)[ind[, 2]]
    
    # Remove dataset.csv file
    file.remove(paste(rootDir, "/Data/Output/dataset.csv", sep = ""))
  }
  
  return(edges)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomiR Number of miRNAs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order miRNAs, mRNAs, TFs
#' @param rootDir Root folder
#' @param usingpVal TRUE if using p-value
#' @param cutoff FDR cutoff
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkWithmiRs = function(interactions, nomiR, nomR, noTF, data, rootDir, usingpVal = FALSE, cutoff = 0.05){
  
  # Build network
  # miRNA => TF
  edges_non_cod <- identifyEdges(data, 1:nomiR, (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir, "miRNATF", usingpVal, cutoff)
  # miRNA => mRNA
  edges_non_cod <- rbind(edges_non_cod,
                         identifyEdges(data, 1:nomiR, (nomiR+1):(nomiR + nomR), rootDir, "miRNAmRNA", usingpVal, cutoff))
  edges_non_cod[, 1] <- gsub("\\.", "-", edges_non_cod[, 1])
  # TF => miRNA
  edges_cod_non <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF), 1:nomiR, rootDir, "TFmiRNA", usingpVal, cutoff)
  edges_cod_non[, 2] <- gsub("\\.", "-", edges_cod_non[, 2])
  # TF => mRNA
  edges_cod_cod <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                                 (nomiR+1):(nomiR + nomR), rootDir, "TFmRNA", usingpVal, cutoff)
  # mRNA => mRNA
  temp <- identifyEdges(data, (nomiR+1):(nomiR + nomR), (nomiR+1):(nomiR + nomR), rootDir, "mRNAmRNA", usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  temp <- identifyEdges(data, (nomiR+nomR+1):(nomiR + nomR + noTF),
                        (nomiR+nomR+1):(nomiR + nomR + noTF), rootDir, "TFTF", usingpVal, cutoff)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_non_cod) <- c("cause", "effect")
  colnames(edges_cod_non) <- c("cause", "effect")
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  # miRNA => TF & miRNA => mRNA: Filter with
  # miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv & TargetScan_7.0.csv
  confirmedList <- read.csv(paste(rootDir, "/Data/miRTarBase_v6.1+TarBase_v7.0+miRWalk_v2.0.csv",
                                  sep = ""), header = FALSE)
  colnames(confirmedList) <- c("cause", "effect")
  predictedList <- read.csv(paste(rootDir, "/Data/TargetScan_7.0.csv", sep = ""), header = TRUE)
  predictedList <- predictedList[, -2]
  colnames(predictedList) <- c("cause", "effect")
  confirmedList[, 1] <- convertmiRs(confirmedList[, 1])
  predictedList[, 1] <- convertmiRs(predictedList[, 1])
  targetbinding <- rbind(confirmedList, predictedList)
  targetbinding <- unique(targetbinding)
  targetbinding$effect <- as.character(targetbinding$effect)
  targetbinding <- targetbinding[which(targetbinding$cause %in% colnames(data)[1:nomiR]),]
  targetbinding <- targetbinding[
    which(targetbinding$effect %in% colnames(data)[(nomiR + 1):(nomiR+nomR+noTF)]),]
  targetbinding <- merge(targetbinding, edges_non_cod)
  # TF => miRNA: Filter with TransmiR from http://www.cuilab.cn/transmir 
  interacts <- read.table(file = paste(rootDir, "/Data/hsa.tsv", sep = ""),
                          sep = '\t', header = FALSE)
  interacts <- interacts[, 1:2]
  colnames(interacts) <- c("cause", "effect")
  interacts[, 2] <- convertmiRs(interacts[, 2])
  interacts <- unique(interacts)
  interacts$cause <- as.character(interacts$cause)
  interacts <- interacts[
    which(interacts$cause %in% colnames(data)[(nomiR+nomR+1):(nomiR+nomR+noTF)]),]
  interacts <- interacts[which(interacts$effect %in% colnames(data)[1:nomiR]),]
  interacts <- merge(interacts, edges_cod_non)
  
  # Combine
  interactions <- rbind(interactions, targetbinding)
  interactions <- rbind(interactions, interacts)
  
  return(interactions)
}

#================================================================
#' Build a network based on the expression data then filtering by existing datasets
#' @param interactions Interactions among mRNAs and TFs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param rootDir Root folder
#' @return Edges of the network from cause to effect
#================================================================
buildNetworkForPersonalised = function(interactions, nomR, noTF, data, rootDir){
  # Build network
  # TF => mRNA
  edges_cod_cod <- identifyEdges(data, (nomR+1):(nomR + noTF),
                                 1:nomR, rootDir)
  # mRNA => mRNA
  temp <- identifyEdges(data, 1:nomR, 1:nomR, rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # TF => TF
  temp <- identifyEdges(data, (nomR+1):(nomR + noTF),
                        (nomR+1):(nomR + noTF), rootDir)
  temp[, 2] <- gsub("\\..*","",temp[, 2])
  edges_cod_cod <- rbind(edges_cod_cod, temp)
  # Set colnames
  colnames(edges_cod_cod) <- c("cause", "effect")
  
  # Filter with existing datasets
  # TF => mRNA, mRNA => mRNA, TF => TF: Filter with interactions
  colnames(interactions) <- c("cause", "effect")
  interactions <- merge(interactions, edges_cod_cod)
  
  return(interactions)
}

#================================================================
#' Analyse a network
#' @param nomiR Number of miRNAs
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param network Edges of the network
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order miRNAs, mRNAs, TFs
#' @param fileName File to be written
#================================================================
analyseNetwork = function(nomiR, nomR, noTF, network, data, fileName){
  # Prepare the data
  dat <- ""
  dat <- paste(dat, "Number of nodes:", nomiR + nomR + noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNAs:", nomiR, "\n", sep = "\t")
  dat <- paste(dat, "Number of TFs:", noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNAs:", nomR, "\n", sep = "\t")
  dat <- paste(dat, "\n", sep = "\t")
  miRs <- colnames(data)[1:nomiR]
  mRs <- colnames(data)[(nomiR+1):(nomiR+nomR)]
  TFs <- colnames(data)[(nomiR+nomR+1):(nomiR+nomR+noTF)]
  edge_type1 <- network[which(network$cause %in% miRs),]
  edge_type1 <- edge_type1[which(edge_type1$effect %in% TFs),]
  edge_type1 <- nrow(edge_type1)
  edge_type2 <- network[which(network$cause %in% miRs),]
  edge_type2 <- edge_type2[which(edge_type2$effect %in% mRs),]
  edge_type2 <- nrow(edge_type2)
  edge_type3 <- network[which(network$cause %in% TFs),]
  edge_type3 <- edge_type3[which(edge_type3$effect %in% miRs),]
  edge_type3 <- nrow(edge_type3)
  edge_type4 <- network[which(network$cause %in% TFs),]
  edge_type4 <- edge_type4[which(edge_type4$effect %in% TFs),]
  edge_type4 <- nrow(edge_type4)
  edge_type5 <- network[which(network$cause %in% TFs),]
  edge_type5 <- edge_type5[which(edge_type5$effect %in% mRs),]
  edge_type5 <- nrow(edge_type5)
  edge_type6 <- network[which(network$cause %in% mRs),]
  edge_type6 <- edge_type6[which(edge_type6$effect %in% mRs),]
  edge_type6 <- nrow(edge_type6)
  edge <- edge_type1 + edge_type2 + edge_type3 + edge_type4 + edge_type5 + edge_type6
  dat <- paste(dat, "Number of edges:", edge, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNA-TF edges:", edge_type1, "\n", sep = "\t")
  dat <- paste(dat, "Number of miRNA-mRNA edges:", edge_type2, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-miRNA edges:", edge_type3, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-TF edges:", edge_type4, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-mRNA edges:", edge_type5, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNA-mRNA edges:", edge_type6, "\n", sep = "\t")
  
  # Write file
  writeLines(dat, fileName)
}

#================================================================
#' Analyse a network
#' @param nomR Number of mRNAs
#' @param noTF Number of transcription factors
#' @param network Edges of the network
#' @param data Expression data with rows being samples and 
#'    columns being biological features with the order mRNAs, TFs
#' @param fileName File to be written
#================================================================
analyseNetworkForPersonalised = function(nomR, noTF, network, data, fileName){
  # Prepare the data
  dat <- ""
  dat <- paste(dat, "Number of nodes:", nomR + noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of TFs:", noTF, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNAs:", nomR, "\n", sep = "\t")
  dat <- paste(dat, "\n", sep = "\t")
  mRs <- colnames(data)[1:nomR]
  TFs <- colnames(data)[(nomR+1):(nomR+noTF)]
  edge_type4 <- network[which(network$cause %in% TFs),]
  edge_type4 <- edge_type4[which(edge_type4$effect %in% TFs),]
  edge_type4 <- nrow(edge_type4)
  edge_type5 <- network[which(network$cause %in% TFs),]
  edge_type5 <- edge_type5[which(edge_type5$effect %in% mRs),]
  edge_type5 <- nrow(edge_type5)
  edge_type6 <- network[which(network$cause %in% mRs),]
  edge_type6 <- edge_type6[which(edge_type6$effect %in% mRs),]
  edge_type6 <- nrow(edge_type6)
  edge <- edge_type4 + edge_type5 + edge_type6
  dat <- paste(dat, "Number of edges:", edge, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-TF edges:", edge_type4, "\n", sep = "\t")
  dat <- paste(dat, "Number of TF-mRNA edges:", edge_type5, "\n", sep = "\t")
  dat <- paste(dat, "Number of mRNA-mRNA edges:", edge_type6, "\n", sep = "\t")
  
  # Write file
  writeLines(dat, fileName)
}

#================================================================
#' Analyse controllability of the network
#' @param result Result file of analysing controllability of the network
#' @param outFile File to be written
#================================================================
analyseControllability = function(result, outFile){
  # Read files
  result_output <- read.csv(result, header = FALSE)
  # Prepare the data
  items <- c("1: # of nodes having links", "2: # of edges", "3: average degree",
             "4: # of driver nodes", "5: fraction of driver nodes = Nd/N",
             "6: fraction of type-I critical nodes",
             "7: fraction of type-I redundant nodes", "8: fraction of type-I ordinary nodes",
             "9: fraction of type-II critical nodes",
             "10: fraction of type-II redundant nodes", "11: fraction of type-II ordinary nodes",
             "12: fraction of critical links",
             "13: fraction of redundant links", "14: fraction of ordinary links")
  dat <- ""
  for (i in 1:14) {
    dat <- paste(dat, items[i], result_output[1,i], "\n", sep = "\t")  
    if(i %in% c(3,5,8,11)) {
      dat <- paste(dat, "\n")
    }
  }
  
  # Write file
  writeLines(dat, outFile)
}

#================================================================
#' Convert a list to a string
#' @param aList A list
#' @return A string contains items of the list
#================================================================
getList = function(aList){
  len <- length(aList)
  l <- aList[1]
  if(len > 1) {
    for(i in 2:len) {
      l <- paste(l, aList[i], sep=",")
    }
  }
  
  return(l)
}

#================================================================
#' Print gene lists for finding overlaps among methods
#' @param OncodriveCLUST A list of genes from OncodriveCLUST method
#' @param OncodriveFM A list of genes from OncodriveFM method
#' @param DawnRank A list of genes from DawnRank method
#' @param CBNA A list of genes from CBNA method
#================================================================
printGeneList = function(OncodriveCLUST, OncodriveFM,
                         DawnRank, CBNA) {
  print("OncodriveCLUST")
  print(getList(OncodriveCLUST))
  print("OncodriveFM")
  print(getList(OncodriveFM))
  print("DawnRank")
  print(getList(DawnRank))
  print("CBNA")
  print(getList(CBNA))  
}

#================================================================
#' Print gene lists for finding overlaps among methods
#' @param OncodriveCLUST A list of genes from OncodriveCLUST method
#' @param ActiveDriver A list of genes from ActiveDriver method
#' @param OncodriveFM A list of genes from OncodriveFM method
#' @param DriverNet A list of genes from DriverNet method
#' @param DawnRank A list of genes from DawnRank method
#' @param CBNA A list of genes from CBNA method
#================================================================
printGeneList_6methods = function(OncodriveCLUST, ActiveDriver, OncodriveFM,
                         DriverNet, DawnRank, CBNA) {
  print("OncodriveCLUST")
  print(getList(OncodriveCLUST))
  print("ActiveDriver")
  print(getList(ActiveDriver))
  print("OncodriveFM")
  print(getList(OncodriveFM))
  print("DriverNet")
  print(getList(DriverNet))
  print("DawnRank")
  print(getList(DawnRank))
  print("CBNA")
  print(getList(CBNA))  
}

#================================================================
#' Print gene lists for finding overlaps among methods
#================================================================
printGeneList_4methods = function(net, DawnRankDriver, NetSig, CBNA) {
  print("DriverNet")
  print(getList(net))
  print("DawnRank")
  print(getList(DawnRankDriver))
  print("NetSig")
  print(getList(NetSig))
  print("CBNA")
  print(getList(CBNA))
}

#================================================================
#' Process data by filtering out genes with constant expression values
#' @param matchedData Matched expression data, rownames are samples and colnames are biological features
#' @return Processed data
#================================================================
processData <- function(matchedData) {
  # Convert input matched expression data to rownames being genes and colnames being samples
  miR <- matchedData$miRs
  mRNA <- matchedData$mRNAs
  Exp = t(cbind(miR, mRNA))
  
  # Filter out genes with constant expression values
  sdGenes <- apply(Exp, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    Exp <- Exp[sdGenes > 0 & !is.na(sdGenes), , drop=FALSE]
  }
  
  return(Exp)
}

#================================================================
#' Classify samples into Epi, Mes, Intermediate Epi, or Intermediate Mes
#' @param Generic_EMT_signature EMT signatures
#' @param Exp Expression data
#' @return Classified samples
#================================================================
scoreEMT <- function(Generic_EMT_signature, Exp) {
  # Get EMT genes, Epi genes, and Mes genes
  EMT_signature = lapply(1:dim(Generic_EMT_signature)[1], function(i) Generic_EMT_signature[i,1])
  Epi = Generic_EMT_signature[which(Generic_EMT_signature[, 2] %in% "Epi"), 1]
  Mes = Generic_EMT_signature[which(Generic_EMT_signature[, 2] %in% "Mes"), 1]
  
  # Only get genes which have expression data
  Epi_Update = Epi[which(Epi %in% rownames(Exp))]
  Mes_Update = Mes[which(Mes %in% rownames(Exp))]
  
  # Estimate GSVA enrichment scores (ES) for two lists of EMT signatures: Epi and Mes.
  # The method is gsva (Hanzelmann et al., 2013), ES is calculated as the maximum distance of the random walk from 0.
  # gsva_es_list = gsva(Exp, list(Epi_Update, Mes_Update), mx.diff=0)$es.obs
  gsva_es_list = gsva(Exp, list(Epi_Update, Mes_Update), mx.diff=0)
  
  # The difference of GSVA enrichment scores between Epi and Mes
  EMT_score = abs(gsva_es_list[2,]) - abs(gsva_es_list[1,])
  
  # Estimate GSVA enrichment score for each gene of EMT signatures.
  # The method is gsva (Hanzelmann et al., 2013), ES is calculated as the maximum distance of the random walk from 0.
  # gsva_es_signature = gsva(Exp, EMT_signature, mx.diff=0)$es.obs
  gsva_es_signature = gsva(Exp, EMT_signature, mx.diff=0)
  
  # Using two-sample KS test, compute the p-values of the difference of GSVA enrichment scores between Epi and Mes
  Epi_index = 1:length(Epi_Update)
  Mes_index = (length(Epi_Update)+1):(length(Epi_Update)+length(Mes_Update))
  EMT_p.value = unlist(lapply(1:dim(gsva_es_signature)[2], function(i) ks.test(gsva_es_signature[Epi_index,i], gsva_es_signature[Mes_index,i])$p.value))
  
  # We divide the samples to four types:
  # Mes (EMT_p.value<0.05 & EMT_score>0),
  # Intermediate Mes (EMT_p.value>=0.05 & EMT_score>0),
  # Epi (EMT_p.value<0.05 & EMT_score<0),
  # and Intermediate Epi (EMT_p.value>=0.05 & EMT_score<0)
  Sample_Phenotype=c()
  for(i in 1:dim(gsva_es_signature)[2]){
    if(EMT_p.value[i]<0.05 & EMT_score[i]>0){
      Sample_Phenotype[i]="Mes"}
    else if(EMT_p.value[i]>=0.05 & EMT_score[i]>0){
      Sample_Phenotype[i]="Intermediate Mes"}
    else if(EMT_p.value[i]<0.05 & EMT_score[i]<0){
      Sample_Phenotype[i]="Epi"}
    else if(EMT_p.value[i]>=0.05 & EMT_score[i]<0){
      Sample_Phenotype[i]="Intermediate Epi"}
  }
  Sample_EMT_Score=cbind(EMT_score,EMT_p.value,Sample_Phenotype)
  
  return(Sample_EMT_Score)
}

#================================================================
#' This function allows you to prepare data for classifying cancer subtypes
#' @param d a matrix of expression of 50 mRNAs in Pam50 with rows being samples and columns being mRNA names.
#' @param directoryPath the directory path to save result file.
#' @param dataDir the directory for saving the result data.
#' @param inputFileName the file name of the input data which is used to classify cancer subtypes.
#================================================================
prepareData=function(d, directoryPath, dataDir = "bioclassifier_data", inputFileName = "inputFile.txt"){
  d <- t(d)
  noOfRow <- nrow(d)
  noOfCol <- ncol(d)
  temp <- matrix(0, nrow = (noOfRow + 1), ncol = noOfCol)
  row.names(temp) <- c("T", row.names(d))
  colnames(temp) <- colnames(d)
  temp[2:(noOfRow + 1), 1:noOfCol] <- d[,]
  write.table(temp, paste(directoryPath, "/", dataDir, "/", inputFileName, sep = ""), col.names=NA, sep = "\t")
}

#================================================================
#' This function allows you to compute measures (Precision, Recall, and F1 Score)
#' @param noOfConfirmed the number of found genes in CGC.
#' @param noOfNotConfirmed the number of found genes not in CGC.
#' @param noOfCGC the number of genes in CGC.
#' @return A list including Precision, Recall, and F1 Score
#================================================================
computeMeasures=function(noOfConfirmed, noOfNotConfirmed, noOfCGC){
  p <- noOfConfirmed/(noOfConfirmed + noOfNotConfirmed)
  r <- noOfConfirmed/noOfCGC
  f1 <- 2*((p*r)/(p+r))
  
  l <- list(Precision = p, Recall = r, F1Score = f1)
  
  return(l)
}

#================================================================
#' This function allows you to prepare data for comparing methods
#' @param Found_Gene_List_Validated the validated found gene list.
#' @param Found_Gene_List the found gene list.
#' @param gold_standard the CGC list.
#' @return A list of genes with corresponding Precision, Recal, and F1 Score for top found genes
#================================================================
prepareDataForComparison=function(Found_Gene_List_Validated, Found_Gene_List, gold_standard){
  noOfCGC <- length(gold_standard)
  noOfFoundGenes <- length(Found_Gene_List)
  r <- matrix(nrow = noOfFoundGenes, ncol = 5)
  colnames(r) <- c("Gene", "InCGC", "Precision", "Recall", "F1Score")
  r[,1] <- Found_Gene_List
  r[,2] <- "No"
  r[which(r[,1] %in% Found_Gene_List_Validated),2] <- "Yes"
  for (i in 1:noOfFoundGenes) {
    topGenes <- r[1:i, , drop=FALSE]
    noOfConfirmed <- sum(topGenes[,2] == "Yes")
    noOfNotConfirmed <- i - noOfConfirmed
    t <- computeMeasures(noOfConfirmed, noOfNotConfirmed, noOfCGC)
    r[i,3] <- t$Precision
    r[i,4] <- t$Recall
    r[i,5] <- t$F1Score
  }
  
  return(r)
}

#================================================================
#' This function allows you to draw a line chart
#' @param OncodriveCLUSTData the data for OncodriveCLUST.
#' @param OncodriveFMData the data for OncodriveFM.
#' @param DawnRankData the data for DawnRank.
#' @param CBNAData the data for CBNA.
#' @param Type it can be "Precision", "Recall", or "F1Score.
#================================================================
drawLineChart=function(OncodriveCLUSTData, OncodriveFMData, DawnRankData, CBNAData, Type="Precision"){
  # Parameters for each chart type
  if(Type == "Precision") {
    yText <- "Precision according to CGC"
    t <- "Precision Comparison"
    pos <- "topright"
  } else if (Type == "Recall") {
    yText <- "Recall according to CGC"
    t <- "Recall Comparison"
    pos <- "topleft"
  } else { # "F1Score"
    yText <- "F1 Score according to CGC"
    t <- "F1 Score Comparison"
    pos <- "topleft"
  }
  
  # Prepare data
  Method=c(rep("OncodriveCLUST", 200), rep("OncodriveFM", 200), rep("DawnRank", 200), rep("CBNA", 200))
  TopNGenes=c(1:200, 1:200, 1:200, 1:200)
  Val=c(OncodriveCLUSTData, OncodriveFMData, DawnRankData, CBNAData)
  indexes<-seq(5,800,5)
  Method <- Method[indexes]
  TopNGenes <- TopNGenes[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
            TopNGenes=TopNGenes,
            Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("OncodriveCLUST", "OncodriveFM", "DawnRank", "CBNA")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$TopNGenes) 
  yrange <- range(data$Val) 
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Top N genes",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], TopNGenes=data$TopNGenes[ind], Val=data$Val[ind])
    lines(method$TopNGenes, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
         pch=plotchar, lty=linetype, title="Method")
}

#================================================================
#' This function allows you to draw a line chart
#' @param OncodriveCLUSTData the data for OncodriveCLUST.
#' @param ActiveDriverData the data for ActiveDriver.
#' @param OncodriveFMData the data for OncodriveFM.
#' @param DriverNetData the data for DriverNet.
#' @param DawnRankData the data for DawnRank.
#' @param CBNAData the data for CBNA.
#' @param Type it can be "Precision", "Recall", or "F1Score.
#================================================================
drawLineChart6=function(OncodriveCLUSTData, ActiveDriverData, OncodriveFMData, DriverNetData, DawnRankData, CBNAData, Type="Precision"){
  # Parameters for each chart type
  if(Type == "Precision") {
    yText <- "Precision according to CGC"
    t <- "Precision Comparison"
    pos <- "topright"
  } else if (Type == "Recall") {
    yText <- "Recall according to CGC"
    t <- "Recall Comparison"
    pos <- "topleft"
  } else { # "F1Score"
    yText <- "F1 Score according to CGC"
    t <- "F1 Score Comparison"
    pos <- "topleft"
  }
  
  # Prepare data
  Method=c(rep("OncodriveCLUST", 200), rep("ActiveDriver", 200), rep("OncodriveFM", 200), rep("DriverNet", 200), rep("DawnRank", 200), rep("CBNA", 200))
  TopNGenes=c(1:200, 1:200, 1:200, 1:200, 1:200, 1:200)
  Val=c(OncodriveCLUSTData, ActiveDriverData, OncodriveFMData, DriverNetData, DawnRankData, CBNAData)
  indexes<-seq(5,1200,5)
  Method <- Method[indexes]
  TopNGenes <- TopNGenes[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
               TopNGenes=TopNGenes,
               Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("OncodriveCLUST", "ActiveDriver", "OncodriveFM", "DriverNet", "DawnRank", "CBNA")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$TopNGenes) 
  yrange <- range(data$Val) 
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Top N genes",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], TopNGenes=data$TopNGenes[ind], Val=data$Val[ind])
    lines(method$TopNGenes, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
         pch=plotchar, lty=linetype, title="Method", bty = "n")
}

#================================================================
#' This function allows you to draw a line chart
#================================================================
drawLineChart7=function(OncodriveCLUSTData, ActiveDriverData, OncodriveFMData, DriverNetData, DawnRankData, NetSigData, CBNAData, Type="Precision"){
  # Parameters for each chart type
  if(Type == "Precision") {
    yText <- "Precision according to CGC"
    t <- "Precision Comparison"
    pos <- "topright"
  } else if (Type == "Recall") {
    yText <- "Recall according to CGC"
    t <- "Recall Comparison"
    pos <- "topleft"
  } else { # "F1Score"
    yText <- "F1 Score according to CGC"
    t <- "F1 Score Comparison"
    pos <- "topleft"
  }
  
  # Prepare data
  Method=c(rep("OncodriveCLUST", 200), rep("ActiveDriver", 200), rep("OncodriveFM", 200), rep("DriverNet", 200), rep("DawnRank", 200), rep("NetSig", 200), rep("CBNA", 200))
  TopNGenes=c(1:200, 1:200, 1:200, 1:200, 1:200, 1:200, 1:200)
  Val=c(OncodriveCLUSTData, ActiveDriverData, OncodriveFMData, DriverNetData, DawnRankData, NetSigData, CBNAData)
  indexes<-seq(5,1400,5)
  Method <- Method[indexes]
  TopNGenes <- TopNGenes[indexes]
  Val <- Val[indexes]
  data <- list(Method=Method,
               TopNGenes=TopNGenes,
               Val=Val)
  
  # Convert to numeric for convenience 
  data$Val <- as.numeric(data$Val)
  
  # Prepare some constants
  methodList <- c("OncodriveCLUST", "ActiveDriver", "OncodriveFM", "DriverNet", "DawnRank", "NetSig", "CBNA")
  nMethods <- length(methodList)
  
  # Get the range for the x and y axis 
  xrange <- range(data$TopNGenes) 
  yrange <- range(data$Val) 
  
  # Set up the plot 
  plot(xrange, yrange, type="n", xlab="Top N genes",
       ylab="", cex.lab=2, cex.axis=2)
  title(ylab=yText, line=2.7, cex.lab=2)
  colors <- rainbow(nMethods) 
  linetype <- c(1:nMethods) 
  plotchar <- seq(18, 18+nMethods, 1)
  
  # Add lines 
  for (i in 1:nMethods) {
    ind <- which(data$Method == methodList[i])
    method <- list(Method=data$Method[ind], TopNGenes=data$TopNGenes[ind], Val=data$Val[ind])
    lines(method$TopNGenes, method$Val, type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i], cex=2) 
  } 
  
  # Add a title 
  title(t, cex.main=2.5)
  
  # Add a legend 
  legend(pos, inset = 0.02, legend=methodList, cex=2, col=colors,
         pch=plotchar, lty=linetype, title="Method", bty = "n")
}

#================================================================
#' This function allows you to prepare data for a clustergram
#' @param dat the data to draw.
#' @param termTop the number of top terms to draw.
#' @param geneTop the number of top genes to draw.
#' @return Data for draw clustergram
#================================================================
prepareDataForClustergram=function(dat, termTop = 10, geneTop  = 20){
  d <- dat
  d <- d[1:termTop,]
  d[,1] <- substr(d[,1], nchar(d[,1])-10, nchar(d[,1])-1)
  genes <- NULL
  geneList <- NULL
  for (i in 1:termTop) {
    genes <- c(genes, unlist(strsplit(d[i,4], ";")))
    geneList <- c(geneList, strsplit(d[i,4], ";"))
  }
  genes <- unique(genes)
  nGenes <- length(genes)
  r <- matrix(0, nrow = termTop+1, ncol = nGenes)
  row.names(r) <- c(d[,1], "total")
  colnames(r) <- genes
  for (i in 1:termTop) {
    ind <- which(colnames(r)[] %in% unlist(geneList[i]))
    r[i,ind] <- 1
  }
  for (i in 1:nGenes) {
    r[termTop+1,i] <- sum(r[1:termTop,i])
  }
  r <- r[,order(r[termTop+1,], decreasing = FALSE)]
  r <- r[1:termTop, (nGenes-geneTop + 1):nGenes]
  r <- t(r)
  newData <- cbind("", r)
  newData[,1] <- rownames(r)
  colnames(newData)[1] <- "Gene"
  
  return(newData)
}

#================================================================
#' This function allows you to draw a clustergram
#' @param dat the data to draw.
#' @param t title
#================================================================
drawClustergram=function(dat, t){
  d <- dat
  d.m <- melt(d)
  d.m <- ddply(d.m, .(variable), transform, rescale = rescale(value))
  p <- ggplot(d.m, aes(variable, Gene)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient(low = "#CCCCCC", high = "#E85642")
  base_size <- 15
  p + theme_grey(base_size = base_size) + 
    labs(title= t, x = "Enriched Terms", y = "Predicted Cancer Drivers") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(legend.position = "none", axis.ticks = element_blank(), 
          axis.text.x = element_text(size = base_size*1.5, 
                                     angle = 285, hjust = 0,
                                     colour = "black"),
          axis.text.y = element_text(size = base_size*1.5,
                                     colour = "black"),
          axis.title = element_text(size=base_size*2),
          plot.title = element_text(size=base_size*2.5))
}

#================================================================
#' This function allows you to read results of exploring drivers for cancer subtypes
#' @param resultDir the directory to get data.
#' @param type the cancer subtype.
#' @return Results of cancer subtype drivers
#================================================================
readResult = function(resultDir, type) {
  # Get the folder name
  outDir <- paste(resultDir, "/", type, sep = "")
  
  # Coding with mutations
  coding_mutations = read.csv(
    file = paste(outDir, "/coding_candidate_cancer_drivers_mutations.csv", sep = ""), as.is = TRUE)
  
  # Coding without mutations
  coding_no_mutations = read.csv(
    file = paste(outDir, "/coding_candidate_cancer_drivers_no_mutations.csv", sep = ""), as.is = TRUE)
  
  # Noncoding
  noncoding = read.csv(
    file = paste(outDir, "/noncoding_candidate_cancer_drivers.csv", sep = ""), as.is = TRUE)
  
  # Combine
  all <- rbind(coding_mutations, coding_no_mutations, noncoding)
  l = list(subType = type, coding_mutations = coding_mutations, coding_no_mutations = coding_no_mutations, noncoding = noncoding, all = all)
  
  return(l)
}

#================================================================
#' This function allows you to read results of exploring drivers for cancer subtypes
#================================================================
readResult02 = function(resultDir, type) {
  # Get the folder name
  outDir <- paste(resultDir, "/", type, sep = "")
  
  # Coding with mutations
  coding_mutations = read.csv(
    file = paste(outDir, "/pVal_coding_candidate_cancer_drivers_mutations.csv", sep = ""), as.is = TRUE)
  
  # Coding without mutations
  coding_no_mutations = read.csv(
    file = paste(outDir, "/pVal_coding_candidate_cancer_drivers_no_mutations.csv", sep = ""), as.is = TRUE)
  
  # Noncoding
  noncoding = read.csv(
    file = paste(outDir, "/pVal_noncoding_candidate_cancer_drivers.csv", sep = ""), as.is = TRUE)
  
  # Combine
  all <- rbind(coding_mutations, coding_no_mutations, noncoding)
  l = list(subType = type, coding_mutations = coding_mutations, coding_no_mutations = coding_no_mutations, noncoding = noncoding, all = all)
  
  return(l)
}

#================================================================
#' This function allows you to get drivers which are specific for a subtype
#' @param basalList Basal drivers.
#' @param her2List Her2 drivers.
#' @param lumAList LumA drivers.
#' @param lumBList LumB drivers.
#' @param normalLikeList Normal-like drivers.
#' @return Different drivers of subtypes
#================================================================
findDiff = function(basalList, her2List, lumAList, lumBList, normalLikeList) {
  basalDiff <- setdiff(basalList$Name, c(her2List$Name, lumAList$Name, lumBList$Name, normalLikeList$Name))
  her2Diff <- setdiff(her2List$Name, c(basalList$Name, lumAList$Name, lumBList$Name, normalLikeList$Name))
  lumADiff <- setdiff(lumAList$Name, c(basalList$Name, her2List$Name, lumBList$Name, normalLikeList$Name))
  lumBDiff <- setdiff(lumBList$Name, c(basalList$Name, her2List$Name, lumAList$Name, normalLikeList$Name))
  normalLikeDiff <- setdiff(normalLikeList$Name, c(basalList$Name, her2List$Name, lumAList$Name, lumBList$Name))
  
  l = list(basalDiff = basalDiff, her2Diff = her2Diff, lumADiff = lumADiff, lumBDiff = lumBDiff, normalLikeDiff = normalLikeDiff)
  
  return(l)
}

#================================================================
#' This function allows you to validate with the CGC
#' @param l List of drivers.
#' @param gold_standard The CGC.
#' @return Validated cancer drivers
#================================================================
validateCGC=function(l, gold_standard) {
  Top50 <- l$coding_mutations[1:50,]
  Top100 <- l$coding_mutations[1:100,]
  Top150 <- l$coding_mutations[1:150,]
  Top200 <- l$coding_mutations[1:200,]
  Top50_Validated <- intersect(Top50[, "Name"], gold_standard)
  Top100_Validated <- intersect(Top100[, "Name"], gold_standard)
  Top150_Validated <- intersect(Top150[, "Name"], gold_standard)
  Top200_Validated <- intersect(Top200[, "Name"], gold_standard)
  
  r <- list(Top50_Validated=Top50_Validated, Top100_Validated=Top100_Validated,
            Top150_Validated=Top150_Validated, Top200_Validated=Top200_Validated)
  
  return(r)
}

#================================================================
#' This function allows you to process miRs of Mes
#' Refer to http://www.mirbase.org/help/nomenclature.shtml for more information
#' @param noncoding_gold_standard_Mes List of miRs.
#' @return Processed miRs
#================================================================
processmiR=function(noncoding_gold_standard_Mes) {
  # Put data to returned variable
  r <- noncoding_gold_standard_Mes
  
  # Get miRs which do not start with "miR"
  ind <- which(substr(r[],1,3) != "miR")
  
  # Add prefix
  r[-ind] <- paste("hsa-", r[-ind], sep = "")
  
  # Get version of miRs
  miRs <- r
  version = checkMiRNAVersion(miRs, verbose=FALSE)
  
  # Convert non mature miRNAs' names to mature names
  miRMature = miRNA_PrecursorToMature(miRs, version=version)
  
  # Convert to version 21
  miRNameList_Accession = miRNA_NameToAccession(miRMature[, 2], version=version)
  miRNameList_Converted = miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v21")
  
  # Process miRNAs' names
  # Update the converted items to the miRNA list
  # If a converted item is "NA", get the corresponding value from the original
  miRs <- miRNameList_Converted[, 2]
  naList <- which(is.na(miRNameList_Converted[, 2]))
  for(i in 1:length(naList)){
    k <- naList[i]
    miRs[k] <- miRMature[k, 1]
  }
  
  return(miRs)
}

#================================================================
#' Calulate probability of degrees
#================================================================
calPro=function(g, cum, mod, typ, nodetype) {
  ver <- nodetype[which(nodetype$TypeI == typ),1]
  ver <- as.character(ver)
  
  dd <- degree.distribution(g, cumulative = cum, mode = mod, v = ver)
  probability <- dd[-1]
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  
  d <- degree(g, mode=mod, v = ver)
  degree <- 1:max(d)
  degree = degree[nonzero.position]
  
  r <- list(pro=probability, deg=degree)
  return(r)
}

#================================================================
#' Get drivers for breast cancer
#================================================================
getBreastDrivers=function(gold_standard) {
  
  r <- NULL
  n <- nrow(gold_standard)
  for (i in 1:n) {
    flag <- FALSE
    if(grepl("breast", as.character(gold_standard[i,2]))) {
      flag <- TRUE
    }
    if(grepl("breast", as.character(gold_standard[i,10]))) {
      flag <- TRUE
    }
    if(grepl("breast", as.character(gold_standard[i,11]))) {
      flag <- TRUE
    }
    if(grepl("breast", as.character(gold_standard[i,12]))) {
      flag <- TRUE
    }
    if(flag) {
      r <- rbind(r, gold_standard[i,])
    }
  }
  
  return(r)
}
