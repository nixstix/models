# here we want to transform our 2020+21 TEST data using our 2020_21 TRAIN MOFA model

source('../scripts/libs.R')
source('../scripts/generic_prediction_functions.R')

# load 2020 data
load('../data/split1/2022_TEST.RData')
rownames(test_data$meta) = test_data$meta$Meta.specimen_id

# load MOFA object:
MOFAobject = load_model('../data/MOFA_models/MOFA2_noScale_train2020_21.hdf5')
MOFAobject

assays = c('geneExp', 'ab',  'cytokineL', 'cell_freq')


## transform 2020 data based on 2021 MOFA model
test_data_tform = mofa_transform(test_data = test_data, assays = assays)
test_data_tform[1:30, 1:5]


# format test matrix:
test_data_tform = format_test_matrix(test_data_tform, test_data)

# baseline day0 
test_data.baseline = test_data_tform[test_data_tform$Meta.timepoint == 0 , ]
rownames(test_data.baseline) = test_data.baseline$Meta.subject_id
colnames(test_data.baseline) = paste('Day0.', colnames(test_data.baseline), sep='') 
test_data.baseline[1:21, 1:5]

# Freq Mn at Day 1
p = test_data_tform[test_data_tform$Meta.timepoint == 1, ]
rownames(p) = p$Meta.subject_id
test_data.baseline$Freq.Monocytes = p[rownames(test_data.baseline), 'Freq.Monocytes']

# CCL3 at Day 3
p = test_data_tform[test_data_tform$Meta.timepoint == 3, ]
rownames(p) = p$Meta.subject_id
test_data.baseline$CCL3 = p[rownames(test_data.baseline), 'Expr.ENSG00000277632.1']


# IgG_PT at Day 14
p = test_data_tform[test_data_tform$Meta.timepoint == 14, ]
rownames(p) = p$Meta.subject_id
test_data.baseline$IgG_PT = p[rownames(test_data.baseline), 'IgG_PT']

# T cell pol Day 28
p = test_data_tform[test_data_tform$Meta.timepoint == 30, ]
rownames(p) = p$Meta.subject_id
test_data.baseline$tcellPol = p[rownames(test_data.baseline), 'TcellPol.pol']
test_data.baseline$tcellPol

save(test_data.baseline, file='../data/MOFA_models/test.RData')


dim(test_data.baseline)
test_data.baseline[1:5, 1:5]

