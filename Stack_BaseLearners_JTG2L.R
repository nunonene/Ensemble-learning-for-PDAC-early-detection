# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

# Stack base-learners

# @author: Dr Nuno R. Nené (Email: nunonene@gmail.com)

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# This code helps replicating the main results of the paper available at https://www.medrxiv.org/content/10.1101/2021.12.02.21267187v2 

# It stacks the base-learners developed with joined time-groups displayed in Figure 1. See Methods section for details.

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

#    In order to reproduce the results displayed in the main text figures data requestors will need to sign a data access agreement and in keeping with patient consent for secondary use obtain ethical approval for any new analyses. 
#    For further inquiries: Dr Nuno R. Nené (Email: nuno.nene.10@ucl.ac.uk)

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Load packages
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

rm(list=ls())

library(caret)
library(pROC)
library(tidyverse)
library(BMA)



# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Helper functions
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

# Extract the prediction probabilities of being a Case from each base-learner model

bestPreds <- function(x){
  stopifnot(is(x, "train"))
  stopifnot({
    x$control$savePredictions %in% c("all", "final") |
      x$control$savePredictions == TRUE
  })
  a <- data.table(x$bestTune, key=names(x$bestTune))
  b <- data.table(x$pred, key=names(x$bestTune))
  b <- b[a, ]
  setorderv(b, c("Resample", "rowIndex"))
  return(b)
}


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Select a particular time-group of interest
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# 2:0-1, 3:0-2, 4:0-3, 5:0-4, 6:0-4+ if focusing on joined time-groups

tg<-2 # Select time-group 


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Load results and stack base-learners
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

BaseLearnerList<-c("C5.0","glmnet","RRF","AdaBag","xgbDART","glmStepAIC","gaussprRadial","svmRadial","pcaNNet","nb")


PredBaseLearners<-lapply(BaseLearnerList,function(algo){
  
  Filename<-paste0('./Results_WithCVandUpSampling_JoinedTimeGroup_',tg,'_',algo,'.Rdata') # Change the file names to match the resampling strategy used for each base learner
  
  if (file.exists(FileName)){load(FileName)}
 
  mat<-bestPreds(trainAlgo) # Extract probability predictions for the best parameters
  
  return(list(pred=mat$Case))# Given that the CV folds were imposed to be the same, the pred should be on the same individuals across algorithms
})


ProbMAT<-data.frame(do.call(cbind,PredBaseLearners)) # Concatenate probability vectors. See schematic in Supplementary Information for illustration purposes

colnames(probMAT)<-BaseLearnerList

# Extract all obs 

load(paste0('./Results_WithUpSampling_',tg,'_TimeGroup',BaseLearnerList[1],'.Rdata')) # Change tg to the desired time-group

mat<-bestPreds(trainAlgo) # Extract probability predictions for the best parameters

status<-factor(as.character(mat$obs),levels=c('Control','Case')) # Given that the CV folds were imposed to be the same, the obs, i.e. real categories, should be the same across algorithms

BMAStack<- bic.glm(x=ProbMAT, y=status, strict = FALSE, OR = 20,glm.family="binomial", factor.type=TRUE) # Fit Bayesian Model Averaging logistic regression model from BMA R package. See paper for package version.


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Load test set
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


TestSet<-read.csv('./TestSet.csv') # The division into Training and Test set was performed by stratifying for a number of clinical covariates and status. See paper for details.
                                    # Data requestors will be provided with the Train-Test set split indices.


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Identify time-groups
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


PCtime_group<-TestSet[,"Time.group"]

UniqueTimeGroups<-as.character(levels(PCtime_group))

TimeGroups_Joined<-lapply(1:(length(UniqueTimeGroups)),function(m){# Joined time-groups=[0-1,0-2,0-3,0-4,0-4+]
  
  if(m==2){
    unlist(lapply(1:m,function(tt){ #0-0.5 and 0.5-1 is taken as one
      which(as.character(PCtime_group)==UniqueTimeGroups[tt])
    }))
  }else{
    
    which(as.character(PCtime_group)==UniqueTimeGroups[m])
  }
  
})

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Create a sub-test for a particular time-group of interest
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# 2:0-1, 3:0-2, 4:0-3, 5:0-4, 6:0-4+ if focusing on joined time-groups

tg<-2 # Select time-group in the test set

PCdataTG<-TestSet[c(TimeGroups_Ind[[tg]]),]



# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Predict probability of being a case for each sample in the test set with the models computed with Run_BaseLearner.R 
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


PredBaseLearners_Test<-lapply(BaseLearnerList,function(algo){
    
    load(paste0('./Results_WithCVandUpSampling_JoinedTimeGroup_',tg,'_',algo,'.Rdata'))# Change tg to the desired time-group
    
    df_test<-PCdataTG[,colnames(trainAlgo$trainingData)[-c(which(colnames(trainAlgo$trainingData)=='.outcome'))]]
    
    pred <- predict(trainAlgo, newdata=df_test, type="prob")
    
     return(pred$Case)
})


PredTest<-do.call(cbind,PredBaseLearners_Test)

colnames(predT)<-c("C5.0","glmnet","RRF","AdaBag","xgbDART","glmStepAIC","gaussprRadial","svmRadial","pcaNNet","nb")

# BMAStack prediction

StackPred <- predict(BMAStack, newdata=predT, type="prob")

# Calculate performances

pROC_BMAStack<- roc(factor(as.character(PCdataTG$status),levels=c('Control','Case')),StackPred)

CI_AUC_TG01_BMA<-ci.auc(pROC_BMAStack,conf.level=0.95, method=c("bootstrap"), boot.n = 2000, boot.stratified = TRUE,
                        progress = getOption("pROCProgress")$name, parallel=TRUE)

Crd<-coords(roc= pROC_BMAStack, x=0.9, input="specificity", transpose = FALSE,ret=c("specificity", "sensitivity","ppv","npv"))

Sensitivity_AtSpe90_BMAStack<-Crd$sensitivity

PPV_AtSpe90_BMAStack<-Crd$ppv 

NPV_AtSpe90_BMAStack<-Crd$npv 


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================
# Save file in order to stack
# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

save(list=c("pROC_BMAStack","Sensitivity_AtSpe90_BMAStack","PPV_AtSpe90_BMAStack","NPV_AtSpe90_BMAStack"),file=paste0('./Results_BMAStack_JTG2L_JoinedTimeGroup_',tg,'.Rdata'))


