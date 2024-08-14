# OforiAnyinam_2024

Metabolic Modeling
Code for metabolic modeling analyses in B Ofori-Anyinam, 2024

Jason H. Yang Lab @ Rutgers New Jersey Medical School

Author: Gautam Mereddy

These scripts perform metabolic modeling analyses on RNA sequencing profiles generated in B Ofori-Anyinam, 2024.

Modeling analyses consisted of two steps:
   1. Creating condition-specific models
   2. Creating flux samples

CREATING CONDITION-SPECIFIC MODELS
Mtb metabolism was modeled using the iEK1011 genome-scale metabolic model (Kavvas ES, BMC System Biology 2018). The De Jesus-Essen version was selected as the base model condition.

RNA sequencing expression profiles were averaged for each biological condition (WT, WT + BDQ, ∆katG, ∆katG + BDQ):
- qs-sum-wt.txt
- qs-sum-bdq.txt
- qs-sum-katg.txt
- qs-sum-katg-bdq.txt

Genes were assigned into high (1), medium (0), or low (-1) expression for each biological condition:
- gd-wt.txt
- gd-bdq.txt
- gd-katg.txt
- gd-katg-bdq.txt

iMAT was applied to generate condition-specific models
- iterate_imat.R
- wt.mat
- bdq.mat
- katg.mat
- katg_bdq.mat

FLUX SAMPLING
Flux variability analysis was performed in CobraPy by using optGpSampler. 10,000 flux samples were collected per condition.
- iterate_pysampling.py
- flux-samples-summary.xlsx
