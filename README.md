# Ensemble-learning-for-PDAC-early-detection


This code helps replicating the main results presented in Figures 1, 2 and 3 of the paper available at https://www.medrxiv.org/content/10.1101/2021.12.02.21267187v2 

It helps to create each base-learner and the stacking strategy with joined time-groups (JTG2L), either with over-sampling of the minority class, synthetic data generation with a High-dimensional Dirichlet Process Mixture (HD-DPM) model or with the Synthetic Minority Oversampling TEchnique (SMOTE):

-- Run_BaseLearners_JoinedTimeGroups.R

-- Run_BaseLearners_WithHDDPM_JoinedTimeGroups.R (calls RunHDDPM_R.py)

-- Run_BaseLearners_WithSMOTE_JoinedTimeGroups.R

-- Stack_BaseLearners_JTG2L.R

A simple script is also provided for model-agnostic feature importance generation and enrichment analysis:

-- FeatureImportance_EnrichmentAnalysis_JoinedTimeGroups.R

See Methods section in the paper for further details.


# Data Availability

For the full data used in the paper, requestors will need to sign a data access agreement and in keeping with patient consent for secondary use obtain ethical approval for any new analyses. For further inquiries: Dr Nuno R. Nené (Email: nuno.nene.10@ucl.ac.uk)
