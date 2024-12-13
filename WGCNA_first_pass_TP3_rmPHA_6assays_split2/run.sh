#!/bin/bash

cd '/home/nthrupp/nthrupp/CMI-PB3/WGCNA_first_pass_TP3_rmPHA_6assays_split2'
pwd

set -e
set -o pipefail

echo 'Creating data directories'

mkdir -p data/MOFA_models
mkdir -p data/split1
mkdir -p data/regression_models

echo 'WGCNA'

cd WGCNA
Rscript first_pass.R &> first_pass.out

cd ../

echo 'Pre-processing'
cd pre-processing


Rscript prep_train_data.R &> prep_train_data.out

cd ../

echo 'MOFA'
cd MOFA

Rscript -e "rmarkdown::render('first_pass_noScale_train.Rmd')"

Rscript prep_testData.R &>  prep_testData.out
Rscript prep_trainData.R &>  prep_trainData.out

cd ../

echo 'Regression models'
cd regression_models

#Rscript predict_CCL3.R &> predict_CCL3.out
Rscript predict_IgG_PT.R &> predict_IgG_PT.out 
Rscript predict_Mn.R &> predict_Mn.out 
Rscript predict_TcellPol.R &> predict_TcellPol.out 
 

#Rscript predict_IgG_PT_1.R &> predict_IgG_PT_1.out 
#Rscript predict_IgG_PT_2.R &> predict_IgG_PT_2.out 

cd ../

echo 'RUN COMPLETE'
