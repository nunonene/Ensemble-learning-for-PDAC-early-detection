#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================

# Run HD-DPM called from R

# @author: Dr Nuno R. Nené (Email: nunonene@gmail.com)

# ===============================================================================================================================================================================
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ===============================================================================================================================================================================


# This code assists in reproducing the main results of the paper available at https://www.medrxiv.org/content/10.1101/2021.12.02.21267187v2 

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



#"""
#Created on Wed Nov 23 16:44:14 2022

#@author: Dr Nuno Rocha Nene
#"""



def Run_DPM_PC(X):

  
    import numpy as np 
    import pandas as pd 
    from sklearn.mixture import BayesianGaussianMixture 


    # Divide into PDAC Cases and Controls

 
    X_Case=X.loc[X['.y'] == 'Case',:]

    X_Control=X.loc[X['.y'] == 'Control',:]
    
    X=X.loc[(X['.y'] == 'Case') | (X['.y'] == 'Control'),:]


    X_Case=X_Case.iloc[:,0:(len(X_Case.iloc[0,])-1)]# Remove status column


    X_Control=X_Control.iloc[:,0:(len(X_Control.iloc[0,])-1)]# Remove status column

    #
    
    X=X.iloc[:,0:(len(X.iloc[0,])-1)]# Remove status column


    #  Fit HD-DPM separately to each class: PDAC and Control


    model_full_Case = BayesianGaussianMixture(
           n_components=len(X_Case.iloc[:,0]), covariance_type='full', weight_concentration_prior=0.01,
           weight_concentration_prior_type='dirichlet_process',
           mean_precision_prior=0.8,
           init_params="kmeans", max_iter=10000, random_state=2) 

    model_full_Control = BayesianGaussianMixture(
           n_components=len(X_Control.iloc[:,0]), covariance_type='full', weight_concentration_prior=0.01,
           weight_concentration_prior_type='dirichlet_process',
           mean_precision_prior=0.8, 
           init_params="kmeans", max_iter=10000, random_state=2) 

    model_full_Case.fit(X_Case);

    model_full_Control.fit(X_Control);
    
 


    max_X_Case=np.zeros(len(X_Case.iloc[0,:]))# Secure that tails of the distribution are trimmed
    min_X_Case=np.zeros(len(X_Case.iloc[0,:]))# Secure that tails of the distribution are trimmed


    for col in range(len(X_Case.iloc[0,:])):
           
        max_X_Case[col]=max(X_Case.iloc[:,col])

        min_X_Case[col]=min(X_Case.iloc[:,col])
    
    
    
    max_X_Control=np.zeros(len(X_Control.iloc[0,:]))# Secure that tails of the distribution are trimmed
    min_X_Control=np.zeros(len(X_Control.iloc[0,:]))# Secure that tails of the distribution are trimmed


    for col in range(len(X_Control.iloc[0,:])):
       
        max_X_Control[col]=max(X_Control.iloc[:,col])

        min_X_Control[col]=min(X_Control.iloc[:,col])
    
   

    max_X=np.zeros(len(X.iloc[0,:]))# Secure that tails of the distribution are trimmed
    min_X=np.zeros(len(X.iloc[0,:]))# Secure that tails of the distribution are trimmed

    for col in range(len(X.iloc[0,:])):
       
        max_X[col]=max(X.iloc[:,col])

        min_X[col]=min(X.iloc[:,col])
    
 
    
    R_ctrl_X=np.ones((100000,len(X_Control.iloc[0,:])))
    R_c_X=np.ones((100000,len(X_Case.iloc[0,:])))


    X_c, y_c = model_full_Case.sample(n_samples=100000)# generate synthetic Cases

    X_ctrl, y_ctrl = model_full_Control.sample(n_samples=100000)# generate synthetic Controls



    for lin in range(len(X_ctrl[:,0])):
    
        for col in range(len(X_ctrl[0,:])):
    
             if X_ctrl[lin,col] < min_X[col]:
                 R_ctrl_X[lin,col]=0
            
             elif X_ctrl[lin,col] > max_X[col]:
         
                 R_ctrl_X[lin,col]=0
            
             if X_c[lin,col] < min_X[col]:
            
                 R_c_X[lin,col]=0
            
             elif X_c[lin,col] > max_X[col]:
         
                 R_c_X[lin,col]=0
                    
            
    y_ctrl=y_ctrl[np.sum(R_ctrl_X,axis=1)==len(X_Control.iloc[0,:])]

    X_ctrl=X_ctrl[np.sum(R_ctrl_X,axis=1)==len(X_Control.iloc[0,:]),:]

    y_c=y_c[np.sum(R_c_X,axis=1)==len(X_Case.iloc[0,:])]

    X_c=X_c[np.sum(R_c_X,axis=1)==len(X_Case.iloc[0,:]),:]


   
    X_HDDPM=pd.DataFrame(np.concatenate((np.hstack((X_c,np.transpose(np.repeat([['Case']], len(X_c[:,0]),axis=1)))),np.hstack((X_ctrl,np.transpose(np.repeat([['Control']], len(X_ctrl[:,0]),axis=1))))),axis=0), columns=np.hstack((np.array(X_Control.columns),np.array('.y'))))
 
    return X_HDDPM


    


