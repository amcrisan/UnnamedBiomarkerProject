library(biom)
library(scales)
library(dplyr)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(e1071)
library(nnet)
library(pROC)
library(vegan)

source("SupportingBiomarkerMethods.R")

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



################################
# Diversity - inverse simpsons
# The AUC and confidence interval measures are slightly off only in the thousandths decimal place
# even though tht exact result is not reproduced I think this is close enough.

# case and non-diarrheal control
tmp <- filter(metadata,disease_stat %in% c("Case","NonDiarrhealControl")) %>%
  select(c(disease_stat,inverseSimpson)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","NonDiarrhealControl"),c(1,0))) %>%
  droplevels

ci.auc(roc(response=tmp$response,predictor=tmp$inverseSimpson)) # 0.808 ( 0.752 - 0.863 ) - agree & slightly off

# case and diarrheal control
tmp <- filter(metadata,disease_stat %in% c("Case","DiarrhealControl")) %>%
  select(c(disease_stat,inverseSimpson)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","DiarrhealControl"),c(1,0))) %>%
  droplevels

tmp$response = factor(tmp$response,levels = c(0,1))
ci.auc(roc(response=tmp$response,predictor=tmp$inverseSimpson)) # 0.581 ( 0.496 - 0.665 ) - agree & slightly - off 


########################################
# The actual biomarker selection
########################################

#grab the LeFSe OTUS from Figure 3
schubertOTUs<-c(1,2,14,22,33,24,10,40,29,12,36,11,15,68,39,
                34,61,38,88,51,63,27,66,46,25,26,37,59,42,56,
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
# some alpha on the Krusal Wallis test - why spend more. The next
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
# penalized regression has also been used in human microarray studies
# and are well suited from when variables >> observations.
#
# Penalized LDA does exist, but I still prefer the more direct interpretations
# of the beta's offorded by glmnet. 

# A. Using the abundance data only
bestOTUs_noMeta<-getBestOTU(response=metadata$disease_stat,
                  countMatrix=abundDat,
                  alph=1,
                  bootNum=100,
                  cvfold=5,
                  logOTUData=TRUE,
                  method = "multinomial")

bestOTUs_noMeta<-melt(table(bestOTUs_noMeta$L1,bestOTUs_noMeta$value))


tmp<-otuCheck(bestOTUs = bestOTUs_noMeta, 
                   taxonomy = taxonomy, 
                   BootMin = 100, 
                   maxTaxaLevel = "all",
                   countMatrix = abundDat,
                   meta  = metadata[,c("sampleID","disease_stat")])

vennList<-vennText(A=schubertOTUs,B=levels(bestOTUs_noMeta$Var.2))



# B. Using metadata
getBestOTU(metadata=metadata.sub,
           response=metadata$disease_stat,
           varsToRemove= NULL,
           countMatrix=abundDat,
           alph=1,
           bootNum=100,
           cvfold=5,
           logOTUData=TRUE,
           method = "multinomial")
