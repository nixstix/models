# REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised by the CMI-PB team. 
# 
# We will take a quick look to see what the data looks like.

## LOAD DATA

ds=c('2020_dataset', '2021_dataset', '2022_dataset', '2023_dataset')

source('../scripts/libs.R')
library(lubridate)

dat = readRDS('../../data/processed_datasets/master_allData_batchCorrected.RDS')
names(dat)

d2 = readRDS('../../data/processed_datasets/master_harmonized_data.RDS')


dob = c(d2$training$subject_specimen$year_of_birth,
  d2$challenge$subject_specimen$year_of_birth)
boost_date = c(d2$training$subject_specimen$date_of_boost, 
        d2$challenge$subject_specimen$date_of_boost)
age = as.numeric(difftime(boost_date,  dob, units = 'weeks')/52)
age

id = c(d2$training$subject_specimen$specimen_id, d2$challenge$subject_specimen$specimen_id)

setdiff(dat$subject_specimen$specimen_id, id)
setdiff(id, dat$subject_specimen$specimen_id)

df = data.frame(specimen_id = id, age = age)
head(df)

dat$subject_specimen = left_join(dat$subject_specimen, df)

## SELECT DS FROM 2020, 2021 AS TRAINING SET, AND 2022 AS TEST SET
dat$subject_specimen = dat$subject_specimen[dat$subject_specimen$dataset %in% ds, ]
table(dat$subject_specimen$dataset)
dim(dat$subject_specimen)

dat$subject_specimen$specimen_id = as.character(dat$subject_specimen$specimen_id)


# MERGE DATA
## merge all matrices to replace empty values with NAs

meta = dat$subject_specimen
colnames(meta) = paste("Meta.", colnames(meta), sep='')
dim(meta)
head(meta)

# ab
x = as.data.frame(t(dat$plasma_ab_titer$batchCorrected_data))
dim(x)
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# cytokines
x = as.data.frame(t(dat$plasma_cytokine_concentrations_by_olink$batchCorrected_data))
colnames(x) = paste("CytokineO.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))


# cytokine legendplex
x = as.data.frame(t(dat$plasma_cytokine_concentrations_by_legendplex$normalized_data))
colnames(x) = paste("CytokineL.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# T cell activation
x = as.data.frame(t(dat$t_cell_activation$raw_data))
colnames(x) = paste("TcellActivation.", colnames(x), sep='')
head(x)
dim(x)
# remove PHA 
x=x[, c(2:3)]
head(x)
dim(x)

meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))



# t cell polarisation
## for now, let's take data from the already-processed data of Pramod
tcellpol = read_tsv('../../data/processed_datasets/Th1(IFN-γ)_Th2(IL-5)_polarization_ratio.tsv')
tcellpol$specimen_id = as.character(tcellpol$specimen_id)
head(tcellpol)

x = tcellpol[,c(1,2)]
head(x)
colnames(x) = c("TcellPol.specimen_id", "TcellPol.pol")
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "TcellPol.specimen_id"))


# cell_freq
x = as.data.frame(t(dat$pbmc_cell_frequency$normalized_data))
colnames(x) = paste("Freq.", colnames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# expression
# WGCNA

load('../data/WGCNA/first_pass_eigengenes-expr.RData')

x = as.data.frame(expr_eigengenes_2$eigengenes)
colnames(x) = paste("Expr.", colnames(x), sep='')
x[1:5, 1:5]
dim(x)
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

dat = list(meta = as.data.frame(meta[,grep('Meta', colnames(meta))]),
                  geneExp = as.data.frame(t(meta[,grep('Expr.', colnames(meta))])) ,
                  ab = as.data.frame(t(meta[,grep('^IgG', colnames(meta))])) ,
                  cytokineO = as.data.frame(t(meta[,grep('CytokineO', colnames(meta))])), 
                  cytokineL = as.data.frame(t(meta[,grep('CytokineL', colnames(meta))])), 
                  tcellPol = as.data.frame(t(meta[,grep('TcellPol', colnames(meta))])),
                  tcellAct = as.data.frame(t(meta[,grep('TcellAct', colnames(meta))])), 
                  cell_freq = as.data.frame(t(meta[,grep('Freq', colnames(meta))])
))

colnames(dat$geneExp) = dat$meta$Meta.specimen_id
colnames(dat$ab) = dat$meta$Meta.specimen_id
colnames(dat$cytokineO) = dat$meta$Meta.specimen_id
colnames(dat$cytokineL) = dat$meta$Meta.specimen_id
colnames(dat$cell_freq) = dat$meta$Meta.specimen_id
colnames(dat$tcellPol) = dat$meta$Meta.specimen_id
colnames(dat$tcellAct) = dat$meta$Meta.specimen_id

dat$ab[1:5, 1:13]

lapply(dat, dim)

# split into test and train

train.idx = dat$meta[dat$meta$Meta.dataset %in% c('2020_dataset', '2021_dataset'), 'Meta.specimen_id']
test.idx = dat$meta[dat$meta$Meta.dataset %in% c('2022_dataset'), 'Meta.specimen_id']
challenge.idx = dat$meta[dat$meta$Meta.dataset %in% c('2023_dataset'), 'Meta.specimen_id']

train_data = list(
  meta = dat$meta[dat$meta$Meta.specimen_id %in% train.idx, ],
  geneExp = dat$geneExp[, train.idx], 
  ab = dat$ab[, train.idx], 
  cytokineO = dat$cytokineO[, train.idx],
  cytokineL = dat$cytokineL[, train.idx],
  cell_freq = dat$cell_freq[, train.idx],
  tcellPol =dat$tcellPol[, train.idx],
  tcellAct = dat$tcellAct[, train.idx]
  
)
train_data$cytokineO
train_data$tcellPol
lapply(train_data, dim)
lapply(train_data, class)

test_data = list(
  meta = dat$meta[dat$meta$Meta.specimen_id %in% test.idx, ],
  geneExp = dat$geneExp[, test.idx], 
  ab = dat$ab[, test.idx], 
  cytokineO = dat$cytokineO[, test.idx],
  cytokineL = dat$cytokineL[, test.idx],
  cell_freq = dat$cell_freq[, test.idx],
  tcellPol =dat$tcellPol[, test.idx],
  tcellAct = dat$tcellAct[, test.idx]
)
test_data$ab[1:5, 1:13]
test_data$tcellPol
lapply(test_data, dim)

challenge_data = list(
  meta = dat$meta[dat$meta$Meta.specimen_id %in% challenge.idx, ],
  geneExp = dat$geneExp[, challenge.idx], 
  ab = dat$ab[, challenge.idx], 
  cytokineO = dat$cytokineO[, challenge.idx],
  cytokineL = dat$cytokineL[, challenge.idx],
  cell_freq = dat$cell_freq[, challenge.idx],
  tcellPol =dat$tcellPol[, challenge.idx],
  tcellAct = dat$tcellAct[, challenge.idx]
)
challenge_data$ab[1:5, 1:13]
challenge_data$tcellPol
lapply(challenge_data, dim)



save(train_data, file='../data/split1/2020+21_TRAIN.RData')
save(test_data, file='../data/split1/2022_TEST.RData')
save(challenge_data, file='../data/split1/2023_CHALLENGE.RData')

