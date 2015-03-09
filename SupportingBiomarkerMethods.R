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




getBestOTU<-function(metadata=NULL,response=NULL,varsToRemove= NULL,countMatrix=NULL,
                     alph=1,bootNum=100,cvfold=4,logOTUData=TRUE,
                     responseType = "binary", ...){
  
  #here is a beautiful vignette on this method:
  #http://web.stanford.edu/~hastie/glmnet/glmnet_beta.html#log
  
  #remove any metadata variables that don't wish to be included
  metadata<-metadata[,setdiff(colnames(metadata),varsToRemove)]
  
  #log the countMatrix data 
  if(logOTUData){
    countMatrix  <- log2(countMatrix + 1)
  }
  
  #adding metadata
  if(is.null(metadata)){
    tmp2 <- countMatrix
  }else{
    #join the metadata with the adbundance data, then drop the sample number
    # TO DO: eventually this has to generalize and maybe not rely on sampleID so much.
    tmp <-merge(x=countMatrix,y=metadata,by.x=0,by.y="sampleID")
    tmp2<-tmp[,2:ncol(tmp)]
  }
  
  
  response <- response
  predictors <- data.matrix(tmp2)
  
  # this is a sanity check that was used for testing
  # which(apply(tmp2,2,function(x){sum(is.na(x))})>0)
  
  
  # run the glment using cv.glment to select the best lambda 
  # default error for logistic regression (multi and bin) : class error
  # default error for all others : deviance (default)
  # need to run cv.glment several times, as it has been documented that on small sample sets
  # and in certain circumstances, the results may be different.
  
  bestOTUs.all<-c()
    
  for(i in 1:bootNum){
    
    #Run the lasso different depending on whether we evaluate binary or pairwise diffs
    if(responseType == "continuous"){
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold)
    }else if (responseType == "binary"){
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold,family="binomial", ...)
    }else if (responseType == "multinomial"){
      # might want to use type.multinomial = "grouped", but this can produce errors (failures
      # to converge) - so for now I am not including it because of difficulty handling those errors.
      lassoResults<-cv.glmnet(x=predictors,y=response,alpha=alph,nfolds=cvfold,
                              family="multinomial", ...)
    }
    
    # This plot is for sanity checking
    #plot(lassoResults)
    #print(lassoResults)
    
    results = coef(lassoResults, s = "lambda.min")
    #bestlambda<-lassoResults$lambda.min
    #results<-predict(lassoResults,s=bestlambda,type="coefficients")
    
    #give me the OTUS baby!
    if(responseType == "multinomial"){
      #extract all the none-zeros OTUs and their beta
      tmp <- sapply(results,function(x){
        bestOTUs<-rownames(x)[which(x !=0)]
        bestOTUs<-bestOTUs[bestOTUs !="(Intercept)"]
        unname(cbind(bestOTUs,x[bestOTUs,]))
        }) 
      
      #make it into a more digestable data frame add add results to growing list
      tmp<-ldply(tmp)
      bestOTUs.all<-rbind(bestOTUs.all,tmp)
      
    }else{
      bestOTUs<-rownames(results)[which(results !=0)]
      bestOTUs.all<-c(bestOTUs.all,unname(cbind(bestOTUs,x[bestOTUs,])))
    }    
  }
  
  if(responseType == "multinomial"){
    colnames(bestOTUs.all)<-c("Response","Predictor","Beta")
  }else{
    colnames(bestOTUs.all)<-c("Predictor","Beta")
  }
  
  #make sure beta is numeric and not a factor!
  bestOTUs.all$Beta<-as.numeric(as.character(bestOTUs.all$Beta))
  
  # output all the OTU results (none-aggregated by mean)
  return(bestOTUs.all)  
}



getEffectScore<-function(OTUdat = NULL, bootMin = NULL,response=NULL, extraFilter = TRUE){
  
  if(is.null(bootMin)){
    stop("Please select the minimum number of times a feature must be selected by glmnet [bootMin]")
  }
  
  if(is.null(response)){
    stop("Please select response type supplied to glment [response]")
  }
  

  # filter OTUs that occur fewer than bootMin times (may have been lucky)
  if(response == "multinomial"){
    tmp<-t(table(OTUdat$Response,OTUdat$Predictor))
    keepOTUs <- rownames(tmp)[rowSums(tmp > bootMin) >= 1]
    
    if(length(keepOTUs) == 0){
      stop("No predictors meet boot min requirements")
    }
    
  }else{
    tmp<-table(OTUdat$Predictor)
    keepOTUs <- names(tmp)[tmp > bootMin]
  }

  OTUdat <- filter(OTUdat,Predictor %in% keepOTUs)

  # summarize the effect size by averaging across bootstrapped iterations
  if(response == "multinomial"){
    betaSummary<-aggregate(Beta ~ Predictor + Response, FUN = summary,data = OTUdat)
  }else{
    betaSummary<-aggregate(Beta ~ Predictor, FUN = summary,data = OTUdat)
  }

  # I think this will be helpful. If the Beta.Min or Beta.Max is not the same sign
  # (both positive, or both negative), then there were iterations when this variables
  # was found to be *both* protective and harmful - which is confusing... so drop that
  # from the list!

  if(extraFilter){
    betaSummary<- mutate(betaSummary, noZero = (Beta[,"Min."] >0) == (Beta[,"Max."] > 0)) %>%
      mutate(meanBeta = Beta[,"Mean"])
    print(sprintf("A total of %d features were dropped because their range of values spanned beta of zero",sum(!(betaSummary$noZero))))
    
    #seems that between multiple computers, the filter statement below either
    #worked or didn't (may depend on dplyr version?). Replaced with something 
    #less elegent.
    #betaSummary<- filter(betaSummary,noZero == TRUE)
    betaSummary<-betaSummary[betaSummary$noZero == TRUE,]
  }

  #create a small plot that shows each OTU and it's effect
  if(response == "multinomial"){
    p<-ggplot(betaSummary,aes(x=Response,y=Predictor))
  }else{
    p<-ggplot(betaSummary,aes(x=1,y=Predictor))
  }

  p <- p + geom_tile(aes(fill = meanBeta),colour="black")+
  scale_fill_gradient2(low = "red",mid = "white",high="blue",na.value="lightgrey")+
  ylab("Log odds")+
  theme_bw()

  #return the plot and a summary of the results
  return(list(effectSizePlot = p, effectSizeSummary = betaSummary))
}




taxaCheck  <- function(otuNames = NULL,taxaLevel = NULL,taxData = NULL){
  otuName<-apply(cbind(otuNames,taxaLevel),1,function(x){
    if(is.na(x[2])){return(NA)}
    return(taxData[[x[1]]][as.numeric(x[2])])})
  
  return(otuName)
}



otuCheck<-function(bestOTUs = NULL, taxonomy = NULL, 
                   maxTaxaLevel = "all",countMatrix = NULL,
                   meta  = NULL,showTopX = NULL){
  
  #setup for maxTaxaLevel 
  linnean = 1:6
  names(linnean)  <- c("kingdom","phylum","class","order","family","genus")
  
  #get the taxa level
  bestOTUs$taxaLevel<-sapply(as.character(bestOTUs$Predictor),function(x){
    getTaxDatLowest(otu=x,taxInfo = taxonomy)
  })
  
  #get taxa name at lowest classifable level
  bestOTUs$taxaName<-taxaCheck(otuNames = as.character(bestOTUs$Predictor),
                                      taxaLevel = bestOTUs$taxaLevel,
                                      taxData = taxonomy)
  
  #look taxa at a specific level
  if(maxTaxaLevel !="all" & maxTaxaLevel %in% names(linnean)){
    maxTaxaLevel <- linnean[maxTaxaLevel]
    bestOTUs<-filter(bestOTUs,taxaLevel >= maxTaxaLevel) %>%
      mutate(concatName = paste(Predictor,taxaName,sep=" - "))
  
    }else{
    bestOTUs <- mutate(bestOTUs,concatName =sprintf("%s (%s)",Predictor,taxaName))
  }
  
  #shared OTUs
  sharedTaxa<-ggplot(data=bestOTUs,aes(x=Response,y=concatName))+
    geom_tile(fill= "black",colour="White")+
    theme_bw()+
    ylab("")+
    xlab("")+
    ggtitle("Shared OTUs among groups")
  
  #looking at the abundance
  if(!(is.null(showTopX))){
    tmp<-countMatrix[,levels(bestOTUs$Predictor)]
    topX<-order(colSums(tmp),decreasing=T)[1:showTopX]
    
    countMatrixDF<-melt(log2(tmp[,topX] + 1))
    
  }else{
    countMatrixDF<-melt(log2(countMatrix[,levels(bestOTUs$Predictor)] + 1))
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


