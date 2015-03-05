######################################################## 
# FUNCTIONS
getTaxDatLowest <- function(otu,taxInfo = NULL)
{
  if(!grepl("Otu",otu)){return(NA)}

  #subsets to the taxonomy coloumns
  x <- taxInfo[[otu]]
  x  <- x[grepl("taxonomy",names(x))]
  idx.lowest <- max(which((x != "unclassified")))
  #return(c(idx.lowest,x[idx.lowest]))
  return(idx.lowest)
}


#remove OTUs below a certain threshold
removeLowAbund<-function(countMatrix = NULL,abundThresh=NULL){
  totalAbundance<-rowSums(as.matrix(countMatrix))
  totalReads<-sum(countMatrix)
  
  abundancePercentage<-totalAbundance/totalReads
  countMatrix<-countMatrix[abundancePercentage>abundThresh,]
  
  #returns the transposed count matrix (sample by OTUs - as nature intended)
  return(t(countMatrix))
}



pairwiseDiffs<-function(someData=NULL)
{
  
  if(is.vector(someData))
  {
    temp <- outer(someData,someData,"-")
    temp <- temp[upper.tri(temp)]    
    return(temp)
    
  } else{
    
    tempMat<-c()
    for(otu in colnames(someData)){
      temp <- outer(as.numeric(someData[,otu]),as.numeric(someData[,otu]),"-")
      temp <- temp[(upper.tri(temp))]
      tempMat<-cbind(tempMat,temp)
    }
    colnames(tempMat)<-colnames(someData)
    return(tempMat)
  }
}



getBestOTU<-function(metadata=NULL,response=NULL,varsToRemove= NULL,countMatrix=NULL,
                     alph=1,bootNum=100,cvfold=4,logOTUData=TRUE,method = "binary"){
  
  #remove any metadata variables that don't wish to be included
  metadata<-metadata[,setdiff(colnames(metadata),varsToRemove)]
  
  #log the countMatrix data 
  if(logOTUData){
    countMatrix  <- log2(countMatrix + 1)
  }
  
  if(is.null(metadata)){
    tmp2 <- countMatrix
    
  }else{
  #join the metadata with the adbundance data, then drop the sample number
  # TO DO: eventually this has to generalize and maybe not rely on sampleID so much.
  tmp <-merge(x=countMatrix,y=metadata,by.x=0,by.y="sampleID")
  tmp2<-tmp[,2:ncol(tmp)]
  }
  
  #for this, I am going to get pairiwse differences between counts and metadata variables
  if(method== "pairwise"){
    #get pairwise differences between all the variables
    tmp2<-pairwiseDiffs(tmp)
    idx.response<-colnames(tmp2) == response
    
    response <- tmp2[,idx.response]
    predictors <- tmp2[,!(idx.response)]
  }else{
    #set up predictors and response using binary response
    response <- response
    predictors <- data.matrix(tmp2)
    
  }
  # this is a sanity check that was used for testing
  # which(apply(tmp2,2,function(x){sum(is.na(x))})>0)
  
  
  #run the lasso (depends on alpha actually, I usually set alpha to 1 which is lasso)
  # use 10 fold cross validation to find the best lambda, and then to select the co-efficients
  bestOTUs.all<-c()
  for(i in 1:bootNum){
    
    #Run the lasso different depending on whether we evaluate binary or pairwise diffs
    if(method == "pairwise"){
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold)
    }else if (method == "binary"){
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold,family="binomial")
    }else if (method == "multinomial"){
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold,family="multinomial")
    }
    
    bestlambda<-lassoResults$lambda.min
    results<-predict(lassoResults,s=bestlambda,type="coefficients")
    
    #give me the OTUS baby!
    if(method == "multinomial"){
      tmp <- sapply(results,function(x){
        bestOTUs<-rownames(x)[which(x !=0)]
        bestOTUs[bestOTUs !="(Intercept)"]}) %>% melt
        
      bestOTUs.all<-rbind(bestOTUs.all,tmp)
      }else{
      bestOTUs<-rownames(results)[which(results !=0)]
      bestOTUs.all<-c(bestOTUs.all,bestOTUs[bestOTUs !="(Intercept)"])
    }    
  }
     
  return(bestOTUs.all)  
}



taxaCheck  <- function(otuNames = NULL,taxaLevel = NULL,taxData = NULL){
  otuName<-apply(cbind(otuNames,taxaLevel),1,function(x){
    if(is.na(x[2])){return(NA)}
    return(taxData[[x[1]]][as.numeric(x[2])])})
  
  return(otuName)
}



otuCheck<-function(bestOTUs = NULL, taxonomy = NULL, BootMin = 100, 
                   maxTaxaLevel = "all",countMatrix = NULL,
                   meta  = NULL,showTopX = NULL){
  
  #setup for maxTaxaLevel 
  linnean = 1:6
  names(linnean)  <- c("kingdom","phylum","class","order","family","genus")
  
  
  
  #remove OTUs that pop up inconsistently (presence was dependent of CV splits)
  bestOTUs<-filter(bestOTUs,BootAttempts == BootMin)
  
  #get the taxa level
  bestOTUs$taxaLevel<-sapply(as.character(bestOTUs$OTU),function(x){
    getTaxDatLowest(otu=x,taxInfo = taxonomy)
  })
  
  #get taxa name at lowest classifable level
  bestOTUs$taxaName<-taxaCheck(otuNames = as.character(bestOTUs$OTU),
                                      taxaLevel = bestOTUs$taxaLevel,
                                      taxData = taxonomy)
  
  #look taxa at a specific level
  if(maxTaxaLevel !="all" & maxTaxaLevel %in% names(linnean)){
    maxTaxaLevel <- linnean[maxTaxaLevel]
    bestOTUs<-filter(bestOTUs,taxaLevel >= maxTaxaLevel) %>%
      mutate(concatName = paste(OTU,taxaName,sep=" - "))
  
    }else{
    bestOTUs <- mutate(bestOTUs,concatName =sprintf("%s (%s)",OTU,taxaName))
  }
  
  #shared OTUs
  sharedTaxa<-ggplot(data=bestOTUs,aes(x=Sample,y=concatName))+
    geom_tile(fill= "black",colour="White")+
    theme_bw()+
    ylab("")+
    xlab("")+
    ggtitle("Shared OTUs among groups")
  
  #looking at the abundance
  if(!(is.null(showTopX))){
    tmp<-countMatrix[,levels(bestOTUs_noMeta$OTU)]
    topX<-order(colSums(tmp),decreasing=T)[1:showTopX]
    
    countMatrixDF<-melt(log2(tmp[,topX] + 1))
    
  }else{
    countMatrixDF<-melt(log2(countMatrix[,levels(bestOTUs$OTU)] + 1))
  }
  colnames(countMatrixDF) <- c("sampleID","OTU","Log2Abundance")
  
  #join with the metadata 
  countMatrixDF<-merge(x=countMatrixDF,y=meta,by="sampleID")
  
  
  #plot the abundance, by case status, for each OTU
  abundance<-ggplot(data = countMatrixDF,aes(x=disease_stat, y = Log2Abundance))+
    geom_boxplot()+
    #geom_point(position = position_jitter(width=0.1,height=0.1))+
    theme_bw()+
    facet_wrap(~OTU)+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0))
  
  return(list(sharedTaxaPlot = sharedTaxa, abundancePlot = abundance))
}



vennText<-function(A=NULL,B=NULL){
  both = intersect(A,B)
  A_Only = setdiff(A,B)
  B_Only = setdiff(B,A)
  
  print("#### TEXTUAL VENN ####")
  print(sprintf("Total Overlapping : %d", length(both)))
  print(sprintf("Unique of list A: %d",length(A_Only)))
  print(" --- ")
  print(A_Only)
  print(" --- ")
  print(sprintf("Unique of list B: %d",length(B_Only)))
  print(B_Only)
  
  return(list(both=both,
              A_only = A_Only,
              B_only = B_Only))

}


