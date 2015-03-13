library(biom)
library(scales)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(e1071)
library(nnet)
library(pROC)
library(vegan)
library(glmnet)

source("SupportingBiomarkerMethods.R")

set.seed(1)
########################################
# Loading and wrangling the data
########################################

# loading the OTU data
dat<- read_biom(biom_file = "data/all.final.an.0.03.subsample.0.03.biom")

abundDat <- t(as.matrix(biom_data(dat))) #rows = Samples, columns = OTU
taxonomy <- observation_metadata(dat)


# loading the metadata seperately
# there needs to be some wrangling to account for samples with not metadata
# and metadata with no associated samples. So just reading the metadata in
# rather than dealing with it in the "make.biom" step.

metadata<-read.csv(file="Data/SchubertMetadata.csv",header=T)

#some manual wraning
#sample DA00754 is only represented once as DA00754.2 whereas it is DA00754
#in the metadata, causing some matching issues. The same is true for DA00939
tempRow<-rownames(abundDat)
tempRow[grep("DA00754.2",rownames(abundDat))]<-"DA00754"
tempRow[grep("DA00939.2",rownames(abundDat))]<-"DA00939"

rownames(abundDat)<-tempRow

# ... the following samples are in the abundace table, but not in the metadata
dropSamps <- setdiff(rownames(abundDat),metadata$sampleID)
# for now, drop samples that don't have associated metadata
abundDat <- abundDat[!(rownames(abundDat) %in% dropSamps), ]

# ... the following have metadata but no associated sample
dropSamps <- setdiff(metadata$sampleID,rownames(abundDat))
metadata  <-  filter(metadata,!(sampleID %in% dropSamps))

# make sure the metadata data and abundance table are in the same
# row order. Not necessary, but I just like it this way
abundDat <- abundDat[match(rownames(abundDat),metadata$sampleID),]

if(sum(!(rownames(abundDat) == as.character(metadata$sampleID))) == 0){
  print("Good to Go")
}


######
# The Wrangling 
# removing the mouse data from the metadata and the abundance matricies
dropSamps<- filter(metadata, host_common_name == "human") %>% select(sampleID)
metadata  <-  filter(metadata,sampleID %in% dropSamps$sampleID) %>% droplevels

abundDat <- abundDat[match(rownames(abundDat),metadata$sampleID),]

########################################
# Sanity check - comparing demographic data to Schubert's paper
########################################
#Initial statistical analyses were conducted to assess differences among 
#the three study groups (C. difficile cases, diarrheal con- trols, and 
#nondiarrheal controls). For continuous variables (e.g., age and weight), 
#one-way analysis of variance was utilized. For categorical variables,
#Pearson’s chi-square test or Fisher’s exact test was performed when 
#expected cell frequencies were less than or equal to 5.

#1. Distribution of case and controls
table(metadata$disease_stat) # checks out

#2. Age
summary(aov(age ~ disease_stat, data = metadata)) # p = 0.034 (agree)

#3. Race
table(metadata$disease_stat,metadata$race)
#reclassify to black, white and Other/unknown
metadata$race2 <- mapvalues(metadata$race,
          from=levels(metadata$race),
          to=c("Other/unknown","Black","Other/unknown","Other/unknown",
               "Other/unknown","White"))

#using chi.sq I can get 0.712 reported in paper, but fishers provides the more exact p-value
fisher.test(table(metadata$disease_stat,metadata$race2)) # p = 0.752 (off - but agree not sig)
chisq.test(table(metadata$disease_stat,metadata$race2)) # p=0.712 (agree - but maybe not right test choice)

#4. Weight
# raw weight values were not given, so it is not possible to perform an anova 
# (metadat) only contains weight categories.

#5. Drug Use
chisq.test(table(metadata$disease_stat,metadata$antibiotics..3mo)) # p<0.0001 (agree)

#6. Others
chisq.test(table(metadata$disease_stat,metadata$antacid)) # p ~ 0.0001 (agree)
chisq.test(table(metadata$disease_stat,metadata$Surgery6mos)) # p<0.0001 (agree)
fisher.test(table(metadata$disease_stat,metadata$historyCdiff)) # p = 0.793 (agree)
fisher.test(table(metadata$disease_stat,metadata$ResidenceCdiff)) # p =  0.593 (agree)
chisq.test(table(metadata$disease_stat,metadata$Healthworker)) # p<= 0.0007 (off - but agree sig)

########################################
# Sanity check - comparing base model performance
########################################
#add the inverse simpsons diversity biomarkers
metadata$inverseSimpson = diversity(abundDat,index="invsimpson")

################################
# base model

#We used age, gender, race, antibiotic use, antacid use, a vegetarian diet, 
#surgery within the past 6 months, a history of CDI, residence with another person who had CDI, 
#and residence with another person who works in a health care setting as
baseModel  <- multinom(disease_stat ~ age + gender + race2 + antibiotics..3mo + antacid + Surgery6mos + 
         historyCdiff + ResidenceCdiff + Healthworker,data = metadata)

metadata$basePredCase <- predict(baseModel,metadata,type="probs")[,"Case"]

# case and non-diarrheal control
tmp <- filter(metadata,disease_stat %in% c("Case","NonDiarrhealControl")) %>%
  select(c(disease_stat,basePredCase)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","NonDiarrhealControl"),c(1,0))) %>%
  droplevels

tmp$response = factor(tmp$response,levels = c(0,1))

ci.auc(roc(response=tmp$response,predictor=tmp$basePredCase)) # 0.894 (0.852-0.936) - agree & slightly off


# case and diarrheal control
metadata$basePredDcontrol <- predict(baseModel,metadata,type="probs")[,"DiarrhealControl"]

tmp <- filter(metadata,disease_stat %in% c("Case","DiarrhealControl")) %>%
  select(c(disease_stat,basePredDcontrol)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","DiarrhealControl"),c(1,0))) %>%
  droplevels

tmp$response = factor(tmp$response,levels = c(0,1))
ci.auc(roc(response=tmp$response,predictor=tmp$basePredDcontrol)) # 0.566 ( 0.481 - 0.650 ) - quite far off - something weird 



# ...... I think I mis-interpreted "nested logit model" after trying to fit one
# without much success. I then realised, that what is meant is "a series of logistic
# models, where a new variable is successively added" (i.e. the nesting is in the variables
# and not the response).  This appraoch seems to work - even so, I don't know if fitting several
# independent logistics is best. I do like the nnet multinom approach, you get  more of a 
# mixed model, which I think is just a nicer way to handle this data response, maybe we should
# talk about that.

comparisons<-rbind(c("Case","DiarrhealControl"),
                   c("Case","NonDiarrhealControl"),
                   c("DiarrhealControl","NonDiarrhealControl"))


# .. getting a storeing the AUCs
storeAUC<-c()
for(comp in 1:nrow(comparisons)){
  comp<-comparisons[comp,]
  metaSubs <- filter(metadata,disease_stat %in% comp)
  #always make sure to select right reference
  if("Case" %in% comp){
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "Case"))
  }else{
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "DiarrhealControl"))
  }
  
  # .... base Model
  baseModel  <- glm(response ~ age + gender + race2 + antibiotics..3mo + antacid + Surgery6mos + 
                           historyCdiff + ResidenceCdiff + Healthworker,
                    family=binomial(link="logit"),
                    data = metaSubs)
  
  metaSubs$prediction  <- predict(baseModel,metaSubs,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC<-rbind(storeAUC,c(paste(comp,collapse=" - "),"base",round(AUCval[2],3),
                  round(AUCval[1],3),
                  round(AUCval[3],3)))
}




################################
# Diversity - inverse simpsons
# THE AUCs while using the "inverse simpsons" measure matched pretty well with some variablity in the
# of 0.1 to 0.2 off the results reported in the paper. This could be because there are different techniques
# to calculate the area under a curve (different approximations) that probably accounts for the difference.
# that they are very much in the right ball park is good.

# .. getting and storing the AUCs
storeAUC_div<-c()
for(comp in 1:nrow(comparisons)){
  comp<-comparisons[comp,]
  metaSubs <- filter(metadata,disease_stat %in% comp)
  print(comp)
  #always make sure to select right reference
  if("Case" %in% comp){
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "Case"))
  }else{
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "DiarrhealControl"))
  }

  # .... metagenomic model - inverse simpsons
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$inverseSimpson))
  
  storeAUC_div<-rbind(storeAUC_div,c(paste(comp,collapse=" - "),"metagenomic-diversity",round(AUCval[2],3),
                             round(AUCval[1],3),
                             round(AUCval[3],3)))
  
  # .... combined model
  
  combinedModel  <- glm(response ~ age + gender + race2 + antibiotics..3mo + antacid + Surgery6mos + 
                      historyCdiff + ResidenceCdiff + Healthworker + inverseSimpson,
                    family=binomial(link="logit"),
                    data = metaSubs)
  
  metaSubs$prediction  <- predict(combinedModel,metaSubs,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_div<-rbind(storeAUC_div,c(paste(comp,collapse=" - "),"combined (metagenomic = diversity)",round(AUCval[2],3),
                             round(AUCval[1],3),
                             round(AUCval[3],3)))

}




########################################
# The actual biomarker selection
########################################

#grab the LeFSe OTUS from Figure 3
schubertOTUs<-c(1,2,14,22,33,24,10,40,29,13,36,11,15,68,39,
                34,61,38,99,51,63,27,66,46,25,26,37,59,42,56,
                85,45,28,41,30,18,21,6,7,3,5,4) %>% 
                formatC(width=4,flag="0") 
schubertOTUs <-paste0("Otu",schubertOTUs)

#1. Filter 1 :  distribution filter
# remove any OTUs where more than 90% of the data is 0
#
# RATIONALE : using a statistical test with a p-value cutoff
# "spends alpha" - this means that the p-value cutoff for subsequent 
# statistical tests needs tobe more stringent for them to really be 
# useful. Also, you don't want to waste your time performing statistical 
# tests on features  that will be unhelpful because there is no information there.
abFact = 0
percentCounts <- apply(abundDat,2,function(x){
  quantile(x,probs = seq(0.1,0.9,by=0.1))["90%"] > abFact
})


#remove those OTUs
abundDat  <-  abundDat[,percentCounts]


#  --- Begin "modified Lefse portion" ---

# 2. Kruskal Wallis Filter
# I can also adjust the p-values after this, but I've choose not to .. for now..
#
# RATIONALE : used in LeFse already, left it in.

otu.pVal<-apply(abundDat,2,function(x){
  kruskal.test(x,g=metadata$disease_stat)$p.value})

abundDat  <-  abundDat[,otu.pVal <0.05]


# --- Normally LefSe does a Wilcox group wise filter here ----
# I really don't think that filter is necessary. You've already spent
# some alpha on the Kruskal Wallis test - why spend more. The next
# step can handle when variables >> observations, so I am going to
# give the next step as much reasonable data as possible.



#4. Penalized regression
#
# RATIONALE : Why not interpret the betas of a regression? Everyone
# knows how to do this, and does it all the time since regressions
# are such a common instrument. Using the glment package, we can
# use penalized regression to perform feature selection. GLMNET
# allows us to use lasso (alpha = 1), ridge (alpha = 0), and
# so called elastic-net (alpha between 0 and 1). By using either of
# the three methods of penalized regression we can control how many
# correlated features (lasso = least; ridge = most) are selected.
# Penalized regression has also been used in human microarray studies
# and are well suited from when variables >> observations.
#
# Also, in practice it seems that there isn't a huge difference between
# using LDA or vanilla binary regression (some sources include the elements of statistical learning - section 4.4.4: 
#  ... p. 105 - "it is generally felt that logistic regression is a safer, more robust bet than the LDA model,
#  ... relying on fewer assumptions. It is our experience that the models give very similar results,
#  ... even when LDA is used inapproperiately, such as with qualtiative predictors.)
# Thus, the substitution of regression for the LDA isn't supposed to be revolutionary. 
# But it does give more flexibility to work with (in my opinion)
#
# Penalized LDA does exist, but I still prefer the more direct interpretations
# of the beta's afforded by glmnet. 



#  ----- A. Using the abundance data only  & Multinomial Response -----
# could be parallelized to improve efficeny.
# using one core, can take ~ 2 - 5 minutes to run based upon boot num
# personally, I prefer a bootnumb of 100, but I am doing 30 to 
# approximately "LefSe's"
bestOTUs_noMeta<-getBestOTU(response=metadata$disease_stat,
                  countMatrix=abundDat,
                  alph=1,
                  bootNum=30,
                  cvfold=5,
                  logOTUData=TRUE,
                  responseType = "multinomial",
                  type.measure="class")


# effect summary variable has two arguments
# ... effectSummary$effectSizePlot which shows the effect log odds beta
# ... effectSummary$effectSizeSummary which is a table that has each OTU, a summary of the beta (mean, min, max ect.)

effectSummary<-getEffectScore(OTUdat = bestOTUs_noMeta,bootMin=(30*0.9),response="multinomial")
bestOTUs_noMeta_pass<-effectSummary$effectSizeSummary

# otuCheck gives me some diagnostic plots regarding the OTUs that have been selected
# 1. sharedTaxaPlot - this is basically a "venn" diagram, but is more readable. It shows
#                     whether an individual OTU is shared between (multinomial) disease states
# 2. abundacePlot - is a little dense. If you want to make it better, set maxTaxLevel to something
#                   else (for example "genus", or "family"). This produces a boxplot of the abundance for each
#                   of the multinomial controls.

tmp<-otuCheck(bestOTUs = bestOTUs_noMeta_pass, 
                   taxonomy = taxonomy, 
                   maxTaxaLevel = "all",
                   countMatrix = abundDat,
                   meta  = metadata[,c("sampleID","disease_stat")],
                response = "multinomial")

# of note  - there is near perfect overlap if I use *all* otus that are found (i.e. effectSize summary bootMin = 0)
# thus, some differences are due to stringency.
vennList<-vennText(A=schubertOTUs,B=unique(as.character(bestOTUs_noMeta_pass$Predictor)))

# .. getting and storing the AUCs
storeAUC_noMeta<-c()
for(comp in 1:nrow(comparisons)){
  comp<-comparisons[comp,]
  metaSubs <- filter(metadata,disease_stat %in% comp)
  print(comp)
  #always make sure to select right reference
  if("Case" %in% comp){
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "Case"))
  }else{
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "DiarrhealControl"))
  }
  
  # .... metagenomic model - biomarkers
  
  biomarker<-log2(abundDat[metaSubs$sampleID,as.character(unique(bestOTUs_noMeta_pass$Predictor))]+1)
  
  biomarkerMeta<-merge(x=metaSubs,y=biomarker,by.x="sampleID",by.y=0)
  
  fla<-as.formula(paste("response ~ ", paste(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))], collapse="+")))
  
  biomarkerFit<-glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(biomarkerFit,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_noMeta<-rbind(storeAUC_noMeta,c(paste(comp,collapse=" - "),"metagenomic-biomarkers",round(AUCval[2],3),
                             round(AUCval[1],3),
                             round(AUCval[3],3)))
  
  # .... combined model
  fla<-as.formula(paste("response ~ ", paste(c(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))],
                                               "age","gender"," race2","antibiotics..3mo","antacid",
                                               "Surgery6mos","historyCdiff","ResidenceCdiff", "Healthworker"), collapse="+")))
  
  
  combinedModel  <- glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(combinedModel,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_noMeta<-rbind(storeAUC_noMeta,c(paste(comp,collapse=" - "),"combined (metagenomic = biomarkers)",round(AUCval[2],3),
                             round(AUCval[1],3),
                             round(AUCval[3],3)))
  
  
  
  # .... metagenomic model - biomarkers#
  biomarker<-log2(abundDat[metaSubs$sampleID,colnames(abundDat) %in% schubertOTUs])
  
  biomarkerMeta<-merge(x=metaSubs,y=biomarker,by.x="sampleID",by.y=0)
  
  fla<-as.formula(paste("response ~ ", paste(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))], collapse="+")))
  
  biomarkerFit<-glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(biomarkerFit,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_noMeta<-rbind(storeAUC_noMeta,c(paste(comp,collapse=" - "),"metagenomic-biomarkers(schubert)",round(AUCval[2],3),
                                           round(AUCval[1],3),
                                           round(AUCval[3],3)))
  
  # .... combined model
  
  
  fla<-as.formula(paste("response ~ ", paste(c(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))],
                                               "age","gender"," race2","antibiotics..3mo","antacid",
                                               "Surgery6mos","historyCdiff","ResidenceCdiff", "Healthworker"), collapse="+")))
  
  
  combinedModel  <- glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(combinedModel,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_noMeta<-rbind(storeAUC_noMeta,c(paste(comp,collapse=" - "),"combined (metagenomic = biomarkers (schubert))",round(AUCval[2],3),
                                           round(AUCval[1],3),
                                           round(AUCval[3],3)))
  
  
}




#here it is with bootMin = 0 (there is much larger overlap if we further filter for reliable OTUs)
# effectSummary<-getEffectScore(OTUdat = bestOTUs_noMeta,bootMin=0,response="multinomial")
# bestOTUs_noMeta_pass<-effectSummary$effectSizeSummary
# 
# vennList<-vennText(A=schubertOTUs,B=unique(as.character(bestOTUs_noMeta_pass$Predictor)))


# ---- Investigate those unique to schubert

# some OTUs were filtered out at the distribution filter step, so
# remove those, from the vennList.
notInAbund<-setdiff(vennList$A_only,colnames(abundDat))
tmpAbund<-melt(abundDat[,setdiff(vennList$A_only,notInAbund)])
colnames(tmpAbund)<-c("sampleID","OTU","Abundance")
tmpAbund <- merge(tmpAbund,metadata[,c("sampleID","disease_stat")],by="sampleID")

#low median abundance?
ggplot(tmpAbund,aes(x=disease_stat,y=log2(Abundance)))+
  geom_boxplot(notch=1)+
  facet_wrap(~OTU)+
  theme_bw()

# ------ investigate those unique to method
#similarity to schubert
# some OTUs were filtered out at the distribution filter step, so
# remove those, from the vennList.
tmpAbund<-melt(abundDat[,vennList$B_only])
colnames(tmpAbund)<-c("sampleID","OTU","Abundance")
tmpAbund <- merge(tmpAbund,metadata[,c("sampleID","disease_stat")],by="sampleID")

# low abundance, too high a reliable on outliers?
ggplot(data=tmpAbund,aes(x=disease_stat,y=log2(Abundance+1)))+
  geom_boxplot()+
  facet_wrap(~OTU)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1))



# ----- B. Using metadata  & Multinomial Response-----
# now we can "adjust" the biomarkers for the impacts of the different metadata variables.
# Loosely, this should mean that biomarkers which are highly correlated with specific
# metadata variables are less likely to be selected. So the biomarkers that drop
# out are likely to provide "value added" information. This still needs to be
# checked afterwards to make sure that's true. Also - because diversity was
# found to be significant, I would also like to adjust for that. 

#create a subsample of the metadata with variables that I want to "adjust" for. 
# I will adjust for all the elements in the base model
metaVars<-c("sampleID",
            "age",
            "gender",
            "race2",
            "antibiotics..3mo",
            "antacid",
            "Surgery6mos",
            "historyCdiff",
            "ResidenceCdiff",
            "Healthworker")

metadata.sub<- metadata[,metaVars]

#now select OTUs again
bestOTUs_Meta<-getBestOTU(metadata=metadata.sub,
           response=metadata$disease_stat,
           varsToRemove= NULL,
           countMatrix=abundDat,
           alph=1,
           bootNum=30,
           cvfold=5,
           logOTUData=TRUE,
           responseType = "multinomial",
           type.measure="class")

effectSummary<-getEffectScore(OTUdat = bestOTUs_Meta,bootMin=(30*0.9),response="multinomial")
bestOTUs_Meta_pass<-effectSummary$effectSizeSummary


#adding another columns to indicate whether they are OTUs or metadata variables
bestOTUs_Meta_pass$biomarkerType  <- ifelse(grepl("Otu",bestOTUs_Meta_pass$Predictor),"OTU","META")

#now just get at the useful OTUs, see how they do.
bestOTUs_Meta_pass_onlyOTUs<-bestOTUs_Meta_pass[bestOTUs_Meta_pass$biomarkerType == "OTU",] %>% droplevels

#check their usefullness
tmp<-otuCheck(bestOTUs = bestOTUs_Meta_pass_onlyOTUs, 
              taxonomy = taxonomy, 
              maxTaxaLevel = "all",
              countMatrix = abundDat,
              meta  = metadata[,c("sampleID","disease_stat")],
              response = "multinomial")

#similarity to schubert (simliar results to OTUs found when not including metadata)
vennList<-vennText(A=schubertOTUs,B=as.character(unique(bestOTUs_Meta_pass_onlyOTUs$Predictor)))

#similarity to prior list (very similar, there are some differences however )
vennList<-vennText(A=as.character(unique(bestOTUs_noMeta_pass$Predictor)),B=as.character(unique(bestOTUs_Meta_pass_onlyOTUs$Predictor)))


# .. getting and storing the AUCs
storeAUC_Meta<-c()
for(comp in 1:nrow(comparisons)){
  comp<-comparisons[comp,]
  metaSubs <- filter(metadata,disease_stat %in% comp)
  print(comp)
  #always make sure to select right reference
  if("Case" %in% comp){
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "Case"))
  }else{
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "DiarrhealControl"))
  }
    
  # .... metagenomic model - biomarkers
  
  biomarker<-log2(abundDat[metaSubs$sampleID,as.character(unique(bestOTUs_Meta_pass_onlyOTUs$Predictor))]+1)
  #biomarker<-abundDat[metaSubs$sampleID,colnames(abundDat) %in% schubertOTUs]
  
  biomarkerMeta<-merge(x=metaSubs,y=biomarker,by.x="sampleID",by.y=0)
  
  fla<-as.formula(paste("response ~ ", paste(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))], collapse="+")))
  
  biomarkerFit<-glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(biomarkerFit,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_Meta<-rbind(storeAUC_Meta,c(paste(comp,collapse=" - "),"metagenomic-biomarkers (meta)",round(AUCval[2],3),
                                           round(AUCval[1],3),
                                           round(AUCval[3],3)))
  
  # .... combined model
  
  
  fla<-as.formula(paste("response ~ ", paste(c(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))],
                                               "age","gender"," race2","antibiotics..3mo","antacid",
                                               "Surgery6mos","historyCdiff","ResidenceCdiff", "Healthworker"), collapse="+")))
  
  
  combinedModel  <- glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(combinedModel,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_Meta<-rbind(storeAUC_Meta,c(paste(comp,collapse=" - "),"combined (metagenomic = biomarkers (meta))",round(AUCval[2],3),
                                           round(AUCval[1],3),
                                           round(AUCval[3],3)))
  
}



# ----- C. Using binomial response & no  -----
# the result that gives the most "overlap" is using not metadata, so I will see what
# perfomring the glment steps with 

#I think the binary and continous response should work, but I feel like I've still
#maybe missed something when mucking about with the multinomial response.


# ----- D. Continuous response  -----
# initially, I had thought to use C.difficle (OTU 19), but since diversity seems
# to be important, I'll try see if there are biomarkers that could "predict"
# diversity.

metaVars<-c("sampleID",
            "age",
            "gender",
            "race2",
            "antibiotics..3mo",
            "antacid",
            "Surgery6mos",
            "historyCdiff",
            "ResidenceCdiff",
            "Healthworker")

metadata.sub<- metadata[,metaVars]

bestOTUs_Meta_continuous<-getBestOTU(metadata=metadata.sub,
                          response=metadata$inverseSimpson,
                          varsToRemove= NULL,
                          countMatrix=abundDat,
                          alph=1,
                          bootNum=30,
                          cvfold=5,
                          logOTUData=TRUE,
                          responseType = "continuous")

effectSummary<-getEffectScore(OTUdat = bestOTUs_Meta_continuous,bootMin=(30*0.9),response="continuous")
bestOTUs_Meta_continuous_pass<-effectSummary$effectSizeSummary


#adding another columns to indicate whether they are OTUs or metadata variables
bestOTUs_Meta_continuous_pass$biomarkerType  <- ifelse(grepl("Otu",bestOTUs_Meta_continuous_pass$Predictor),"OTU","META")

#a quick look at the useful metadata variables
#filter(bestOTUs_Meta_continuous_pass,biomarkerType == "META")

#now just get at the useful OTUs, see how they do.
bestOTUs_Meta_continuous_pass_onlyOTUs<-bestOTUs_Meta_continuous_pass[bestOTUs_Meta_continuous_pass$biomarkerType == "OTU",] %>% droplevels

tmp<-otuCheck(bestOTUs = bestOTUs_Meta_continuous_pass_onlyOTUs, 
              taxonomy = taxonomy, 
              maxTaxaLevel = "all",
              countMatrix = abundDat,
              meta  = metadata[,c("sampleID","disease_stat")],
              response = "continuous")


# .. getting and storing the AUCs
storeAUC_Meta_div<-c()
for(comp in 1:nrow(comparisons)){
  comp<-comparisons[comp,]
  metaSubs <- filter(metadata,disease_stat %in% comp)
  print(comp)
  #always make sure to select right reference
  if("Case" %in% comp){
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "Case"))
  }else{
    metaSubs<-mutate(metaSubs, response = as.numeric(disease_stat == "DiarrhealControl"))
  }
  
  # .... metagenomic model - biomarkers
  
  biomarker<-log2(abundDat[metaSubs$sampleID,as.character(unique(bestOTUs_Meta_continuous_pass_onlyOTUs$Predictor))]+1)
  #biomarker<-abundDat[metaSubs$sampleID,colnames(abundDat) %in% schubertOTUs]
  
  biomarkerMeta<-merge(x=metaSubs,y=biomarker,by.x="sampleID",by.y=0)
  
  fla<-as.formula(paste("response ~ ", paste(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))], collapse="+")))
  
  biomarkerFit<-glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(biomarkerFit,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_Meta_div<-rbind(storeAUC_Meta_div,c(paste(comp,collapse=" - "),"metagenomic-biomarkers (diversity)",round(AUCval[2],3),
                                       round(AUCval[1],3),
                                       round(AUCval[3],3)))
  
  # .... combined model
  
  
  fla<-as.formula(paste("response ~ ", paste(c(colnames(biomarkerMeta)[grepl("Otu",colnames(biomarkerMeta))],
                                               "age","gender"," race2","antibiotics..3mo","antacid",
                                               "Surgery6mos","historyCdiff","ResidenceCdiff", "Healthworker"), collapse="+")))
  
  
  combinedModel  <- glm(fla,data=biomarkerMeta,family=binomial(link="logit"),maxit=1000)
  
  metaSubs$prediction  <- predict(combinedModel,biomarkerMeta,type="response")
  
  AUCval<-ci.auc(roc(response=metaSubs$response,predictor=metaSubs$prediction))
  
  storeAUC_Meta_div<-rbind(storeAUC_Meta_div,c(paste(comp,collapse=" - "),"combined (metagenomic = biomarkers (diversity))",round(AUCval[2],3),
                                       round(AUCval[1],3),
                                       round(AUCval[3],3)))
  
}





######################################
# combine all of the AUC tables floating around

tab<-rbind(storeAUC,
storeAUC_div,
storeAUC_noMeta,
storeAUC_Meta,
storeAUC_Meta_div)

write.csv(tab,file="Reports/AUC_for_metadata.csv",quote=F)
