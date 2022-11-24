# Ensemble-learning-for-PDAC-early-detection


# This code helps replicating the main results of the paper available at https://www.medrxiv.org/content/10.1101/2021.12.02.21267187v2 

# It helps generating each base-learner and the stacking strategy with joined time-groups JTG2L with over-sampling of the minority class and synthetic data generation with HD-DPM or SMOTE are available 

# A simple script is also provided for model-agnostic feature importance generation and enrichment analysis


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

#    For full the data used in the paper, requestors will need to sign a data access agreement and in keeping with patient consent for secondary use obtain ethical approval for any new analyses. 
#    For further inquiries: Dr Nuno R. Nené (Email: nuno.nene.10@ucl.ac.uk)
