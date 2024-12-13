

source('../scripts/libs.R')
source('../scripts/generic_prediction_functions.R')

MOFAobject = load_model('../data/MOFA_models/MOFA2_noScale_train2020_21.hdf5')
MOFAobject

# get predictors
# --------------

predictors = format_training_matrix(mofa_object = MOFAobject) # concatenate MOFA factors, clinical data and assay data ; code categorical variables to numeric
predictors[1:5, 1:5]
predictors[1:5, 5110:ncol(predictors)]

# extract day 0 
# --------------

predictors.baseline = predictors[predictors$Meta.timepoint == 0, ]
rownames(predictors.baseline) = predictors.baseline$Meta.subject_id
colnames(predictors.baseline) = paste('Day0.', colnames(predictors.baseline), sep='') 
dim(predictors.baseline)
predictors.baseline[1:5, 1:5]

# extract tasks
# -------------

load('../data/split1/2020+21_TRAIN.RData')


# Freq Mn at Day 1
p = train_data$meta[train_data$meta$Meta.timepoint == 1, ]
p$Freq.Monocytes = as.numeric(train_data$cell_freq['Freq.Monocytes',p$Meta.specimen_id])
rownames(p) = p$Meta.subject_id
predictors.baseline$Freq.Monocytes = p[rownames(predictors.baseline), 'Freq.Monocytes']

# CCL3 at Day 3
p = train_data$meta[train_data$meta$Meta.timepoint == 3, ]
p$Expr.ENSG00000277632.1 = as.numeric(train_data$geneExp['Expr.ENSG00000277632.1',p$Meta.specimen_id]) 
rownames(p) = p$Meta.subject_id
predictors.baseline$CCL3 = p[rownames(predictors.baseline), 'Expr.ENSG00000277632.1']


# IgG_PT at Day 14
p = train_data$meta[train_data$meta$Meta.timepoint == 14, ]
p$IgG_PT = as.numeric(train_data$ab['IgG_PT',p$Meta.specimen_id]) 
rownames(p) = p$Meta.subject_id
predictors.baseline$IgG_PT = p[rownames(predictors.baseline), 'IgG_PT']

# T cell activation Day 28
p = train_data$meta[train_data$meta$Meta.timepoint == 30, ]
p$tcellPol = as.numeric(train_data$tcellPol['TcellPol.pol',p$Meta.specimen_id]) 
rownames(p) = p$Meta.subject_id
predictors.baseline$tcellPol = p[rownames(predictors.baseline), 'tcellPol']
predictors.baseline$tcellPol

# save
# ----

save(predictors.baseline, file='../data/MOFA_models/predictors.RData')

