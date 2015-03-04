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

# ... the following samples are in the abundace table, but not in the metadata
dropSamps <- setdiff(rownames(abundDat),metadata$sampleID)
# for now, drop samples that don't have associated metadata
abundDat <- abundDat[!(rownames(abundDat) %in% dropSamps), ]

# ... the following have metadata but no associated sample
dropSamps <- setdiff(metadata$sampleID,rownames(abundDat))
metadata  <-  filter(metadata,!(sampleID %in% dropSamps))

######
# The Wrangling 
# removing the mouse data from the metadata and the abundance matricies
dropSamps<- filter(metadata, host_common_name == "human") %>% select(sampleID)
metadata  <-  filter(metadata,sampleID %in% dropSamps$sampleID) %>% droplevels

abundDat <- abundDat[match(rownames(abundDat),metadata$sampleID),]

#sanity check - make sure the data are in teh same older
sum(rownames(abundDat) != metadata$sampleID) == 0


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
table(metadata$disease_stat) # I have 2 fewer cases

#2. Age

#3. Race
table(metadata$disease_stat,metadata$race)
#reclassify to black, white and Other/unknown
metadata$race2 <- mapvalues(metadata$race,
          from=levels(metadata$race),
          to=c("Other/unknown","Black","Other/unknown","Other/unknown",
               "Other/unknown","White"))

fisher.test(table(metadata$disease_stat,metadata$race2)) # p = 0.79 (off - but agree not sig)

#4. Weight

#5. Drug Use
chisq.test(table(metadata$disease_stat,metadata$antibiotics..3mo)) # p<0.0001 (agree)

#6. Others
chisq.test(table(metadata$disease_stat,metadata$antacid)) # p ~ 0.0001 (agree)
chisq.test(table(metadata$disease_stat,metadata$Surgery6mos)) # p<0.0001 (agree)
fisher.test(table(metadata$disease_stat,metadata$historyCdiff)) # p = 0.793 (agree)
fisher.test(table(metadata$disease_stat,metadata$ResidenceCdiff)) # p =  0.595 (off - but agree not sig)
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

ci.auc(roc(response=tmp$response,predictor=tmp$basePredCase)) # 0.891 (0.849-0.934) - agree & slightly off


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

# case and non-diarrheal control
tmp <- filter(metadata,disease_stat %in% c("Case","NonDiarrhealControl")) %>%
  select(c(disease_stat,inverseSimpson)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","NonDiarrhealControl"),c(1,0))) %>%
  droplevels

ci.auc(roc(response=tmp$response,predictor=tmp$inverseSimpson)) # 0.809 ( 0.754 - 0.866 ) - agree & exact

# case and diarrheal control
tmp <- filter(metadata,disease_stat %in% c("Case","DiarrhealControl")) %>%
  select(c(disease_stat,inverseSimpson)) %>%
  mutate(response = mapvalues(disease_stat,c("Case","DiarrhealControl"),c(1,0))) %>%
  droplevels

tmp$response = factor(tmp$response,levels = c(0,1))
ci.auc(roc(response=tmp$response,predictor=tmp$inverseSimpson)) # 0.583 ( 0.498 - 0.668 ) - agree & slightly - off 


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
abFact = 0
percentCounts <- apply(abundDat,2,function(x){
  quantile(x,probs = seq(0.1,0.9,by=0.1))["90%"] > abFact
})


abundDat  <-  abundDat[,percentCounts]


#  --- Begin "modified Lefse portion" ---

# 2. Kruskal Wallis Filter
# I can also adjust the p-values after this, but I've choose not to .. for now..

otu.pVal<-apply(abundDat,2,function(x){
  kruskal.test(x,g=metadata$disease_stat)$p.value})

abundDat  <-  abundDat[,otu.pVal <0.05]


# --- Normally LefSe does a Wilcox group wise filter here ----
# I really don't think that filter is necessary. You've already spen
# some alpha on the Krusal Wallis test. 



#4. Penalized regression

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
