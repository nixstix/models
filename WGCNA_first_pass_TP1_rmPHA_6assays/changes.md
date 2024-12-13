!!!! still to complete
!!!! prep after MOFA
!!!! select only factors

The original pipeline is specified under `MOFA_7assays`. Any changes to this pipeline are documented below.

## MOFA/first_pass_noScale_train.Rmd

Changed path in first chunk

```
knitr::opts_knit$set(root.dir = '/mnt/BioAdHoc/Groups/Peters/nthrupp/CMI-PB3/MOFA_7assays_TP0/MOFA/')
nr_factors = 15

assays = c("geneExp",  "ab" , "cytokineO" ,"cytokineL", "cell_freq" ,"tcellPol" ,"tcellAct")
```

Subsetted data based on TP0

### subset for only TP0

```
sample_idx = train_data$meta[train_data$meta$Meta.timepoint == 0 , 'Meta.specimen_id'] 
length(sample_idx)

train_data$meta = train_data$meta[train_data$meta$Meta.specimen_id %in% sample_idx, ]
train_data$geneExp = train_data$geneExp[,sample_idx]
train_data$ab = train_data$ab[,sample_idx]
train_data$cytokineO = train_data$cytokineO[,sample_idx]
train_data$cytokineL = train_data$cytokineL[,sample_idx]
train_data$cell_freq = train_data$cell_freq[,sample_idx]
train_data$tcellPol = as.data.frame(train_data$tcellPol[,sample_idx])
train_data$tcellAct = train_data$tcellAct[,sample_idx]
```


## Regression_models/*.R

Changed from model3 (factors+bl) to model2 (factors) as baseline predictive models:

```
# select top features from above
to_test = model2$model_coef
to_test = to_test[grep('Factor', names(to_test))]
to_test = to_test[order(abs(to_test), decreasing = T)]
to_test

```


