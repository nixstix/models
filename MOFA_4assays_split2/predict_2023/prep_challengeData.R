# here we want to transform our 2020+21 TEST data using our 2020_21 TRAIN MOFA model

source('../scripts/libs.R')
source('../scripts/generic_prediction_functions.R')

# create data dir
dir.create('../data/predictions')

# load 2020 data
load('../data/split1/2023_CHALLENGE.RData')
rownames(challenge_data$meta) = challenge_data$meta$Meta.specimen_id

# load MOFA object:
MOFAobject = load_model('../data/MOFA_models/MOFA2_noScale_train2020_21.hdf5')
MOFAobject

assays = c('geneExp', 'ab', 'cytokineL', 'cell_freq')


## transform 2020 data based on 2021 MOFA model
challenge_data_tform = mofa_transform(test_data = challenge_data, assays = assays)
challenge_data_tform[1:30, 1:5]


# format test matrix:
challenge_data_tform = format_test_matrix(challenge_data_tform, challenge_data)

# baseline day0 
challenge_data.baseline = challenge_data_tform[challenge_data_tform$Meta.timepoint == 0 , ]
rownames(challenge_data.baseline) = challenge_data.baseline$Meta.subject_id
colnames(challenge_data.baseline) = paste('Day0.', colnames(challenge_data.baseline), sep='') 
challenge_data.baseline[1:21, 1:5]


save(challenge_data.baseline, file='../data/predictions/challenge.RData')


dim(challenge_data.baseline)
challenge_data.baseline[1:5, 1:5]

