# QUICK LOOK AT DATA
# -------------------------------------------------

# This data has already been normalised and batch-corrected by the CMI-PB team. 
# 

## LOAD DATA

source('../scripts/libs.R')

dat = readRDS('../../data/processed_datasets/master_allData_batchCorrected.RDS')
names(dat)

m = dat$subject_specimen
table(m$dataset, m$subject_id)
# donors 2020 : 60
# donors 2021 : 36
# donors 2022 : 21
# donors 2023 (test) : 55
