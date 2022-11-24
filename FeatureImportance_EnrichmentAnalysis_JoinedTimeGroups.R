rm(list=ls())

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

# Determine feature importance for each base-learner and perform enrichment analysis

# @author: Dr Nuno R. Nené (Email: nunonene@gmail.com)

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# This code helps replicating the main results of the paper available at https://www.medrxiv.org/content/10.1101/2021.12.02.21267187v2 

# It calculates the model agnostic feature importance associated with each the base-learner with joined time-groups, retains only those with non-zero values and performs enrichment analysis with Gprofiler.

# See Methods section and Supplementary Information for details

# Early detection of pancreatic ductal adenocarcinomas with an ensemble learning model based on a panel of protein serum biomarkers

# Authors:  Nuno R. Nené1,*, Alexander Ney2, Tatiana Nazarenko1,3, Oleg Blyuss1,4,5, Harvey E. Johnston1,6, Harry J. Whitwell1,4,7,8, 
#           Eva Sedlak1, Aleksandra Gentry-Maharaj9, Sophia Apostolidou9, Eithne Costello10, William Greenhalf11, Ian Jacobs1,12, Usha Menon9, 
#            Justin Hsuan2, Stephen P. Pereira2, Alexey Zaikin1,3,4,13, John F. Timms1.

# Affiliations:
#   1 Department of Women’s Cancer, EGA Institute for Women’s Health, University College London, 84-86 Chenies Mews, London, WC1E 6HU, UK. 
#   2 Institute for Liver and Digestive Health, University College London, Upper 3rd Floor, Royal Free Campus, Rowland Hill Street, London NW3 2PF, UK.
#   3 Department of Mathematics, University College London, London WC1H 0AY, UK.
#   4 Department of Applied Mathematics, Lobachevsky University of Nyzhniy Novgorod, Nizhniy Novgorod 603105, Russia.
#   5 Wolfson Institute of Population Health, Queen Mary University of London, Charterhouse Square, EC1M 6BQ, London, UK
#   6 Babraham Institute, Babraham Research Campus Cambridge, CB22 3AT, UK
#   7 National Phenome Centre and Imperial Clinical Phenotyping Centre, Department of Metabolism, Digestion and Reproduction, IRDB Building, Imperial College London, Hammersmith Campus, London, W12 0NN, UK
#   8 Section of Bioanalytical Chemistry, Division of Systems Medicine, Department of Metabolism, Digestion and Reproduction, Imperial College London, South Kensington Campus, London, SW7 2AZ, UK
#   9 MRC Clinical Trials Unit at UCL, Institute of Clinical Trials and Methodology, UCL, 90 High Holborn, 2nd Floor, London, WC1V 6LJ, UK
#   10 Department of Molecular and Clinical Cancer Medicine, University of Liverpool, Liverpool, UK
#   11 Liverpool Experimental Cancer Medicine Centre, University of Liverpool, Liverpool, L69 3GL, UK
#   12 University of New South Wales, Sydney, NSW, 2052, Australia
#   13 World-Class Research Center "Digital biodesign and personalized healthcare", Sechenov First Moscow State Medical University, Moscow, Russia 

#   *To whom correspondence should be addressed: 
#    Dr Nuno R. Nené (Email: nuno.nene.10@ucl.ac.uk)


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Data Availability
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

#    Data requestors will need to sign a data access agreement and in keeping with patient consent for secondary use obtain ethical approval for any new analyses. 
#    For further inquiries: Dr Nuno R. Nené (Email: nuno.nene.10@ucl.ac.uk)

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Load packages
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================



library(caret)
library(tidyverse)
library(vip) # For the model agnostic feature importance method
library(gprofiler2) # For enrichment analysis
library(stringr)


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Select a particular time-group of interest
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# 2:0-1, 3:0-2, 4:0-3, 5:0-4, 6:0-4+ if focusing on joined time-groups

tg<-6 # Select time-group. Figure 3 corresponds to finding the enrichment terms overexpressed given the set of features with importance non-zero for the models developed with 0-4+ samples. See Methods section of the paper for further details.


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Load results and calculate feature importance for each base-learner
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

BaseLearnerList<-c("C5.0","glmnet","RRF","AdaBag","xgbDART","glmStepAIC","gaussprRadial","svmRadial","pcaNNet","nb")


FeatImportance<-lapply(BaseLearnerList,function(m){
  
  load(paste0('./Results_WithUpSampling_',tg,'_TimeGroup_',algo,'.Rdata')) # Change tg to the desired time-group. Change file name to that generated with the desired resampling strategy
  
  v<-vi(trainAlgo, method = "firm")

  return(data.frame(Features=v$Variable[which(v$Importance!=0)],Importance=as.matrix(v$Importance[which(v$Importance!=0)]/max(v$Importance[which(v$Importance!=0)])),Algo=rep(m,length(which(v$Importance!=0)))))
  
})


FeatImportance<-do.call(rbind,FeatImportance)

FeatImportance$Features<-as.character(FeatImportance$Features)




# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Save Feature Importance File 
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


save(list=c("FeatImportance"),file=paste0('./FeatureImportance_Results_WithUpSampling_JoinedTimeGroup_',tg,'.Rdata'))# Change filename to match the resampling strategy



# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Enrichment analysis with Gprofiler
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

UniqueFeatures<-unique(as.character(FeatImportance$Features)) # Remove any features that are not proteins and make sure the names are standard by resorting to library(org.Hs.eg.db)


gp<-gost(UniqueFeatures, organism = "hsapiens",evcodes =TRUE)


# Select terms

db<-'KEGG' # e.g.KEGG, REAC, WP, GO:BP,etc

Terms<-which(gp$result$source==db)


TermsPval<-data.frame(Term=str_to_upper(gp$result$term_name[Terms], locale = "en"),Pvalue=-log(gp$result$p_value[Terms],10))



# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Save Feature Importance File 
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


save(list=c("TermsPval"),file=paste0('./FeatureGprofiler_Results_WithUpSampling_JoinedTimeGroup_',tg,'.Rdata'))# Change filename to the desired resampling strategy


