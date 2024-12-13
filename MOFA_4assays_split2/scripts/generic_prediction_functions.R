# function to format features ready for prediction:
# extract MOFA factors from mofa object
# add clinical data and assay data
# recode categorical data to numeric

format_training_matrix = function(mofa_object){
    # factors calculated by MOFA:
    predictors = as.data.frame(MOFAobject@expectations$Z[[1]])
    predictors$Meta.specimen_id = rownames(predictors)
    
    # add assay data
    predictors = merge(predictors, MOFAobject@samples_metadata, by = 'Meta.specimen_id')
    
    rn = predictors$Meta.specimen_id
    
    # categorical to numerical data
    predictors$infancy_vac = recode(predictors$Meta.infancy_vac, 
                                    'aP' = 0, 
                                    'wP' = 1)
    predictors$biological_sex = recode(predictors$Meta.biological_sex, 
                                       "Male" = 0, 
                                       "Female" = 1)
    
    #predictors$Meta.year_of_birth = 2023 - (as.numeric(sapply(strsplit(predictors$Meta.year_of_birth, '-'), `[`, 1)))
    
    rownames(predictors) = rn
    return(predictors)    
}

# function to format features ready for prediction:
# recode categorical data to numeric

format_test_matrix = function(test_data_f, original_data){
    
    test_data_f$Meta.specimen_id = original_data$meta$Meta.specimen_id
    rn = test_data_f$Meta.specimen_id
    
    test_data_f = cbind(test_data_f, t(original_data$geneExp), t(original_data$ab), 
                  t(original_data$cytokineL), t(original_data$cell_freq), t(original_data$cytokineO), 
                  t(original_data$tcellAct), t(original_data$tcellPol), original_data$meta)
    
    test_data_f$infancy_vac = recode(test_data_f$Meta.infancy_vac, 
                               "aP"=0, 
                               "wP"=1)
    test_data_f$Meta.timepoint = test_data_f$Meta.timepoint
    
    
    #x=test_data_f$Meta.year_of_birth
    #x = 2023 - (as.numeric(sapply(strsplit(as.character(x), '-'), `[`, 1)))
    #test_data_f$Meta.year_of_birth = x
    
    rownames(test_data_f) = rn
    return(test_data_f)    
}


# function to transform testing data based on MOFA model
mofa_transform = function(test_data, assays){
  # calculate a weight for each factor in each view / assay
  # we will sum each factor across the 4 views / assays, first weighting each view by its variance explained
  varExpl = calculate_variance_explained(MOFAobject, views = "all", factors = "all")
  varExpl
  weights = t(varExpl$r2_per_factor$group1) / 100
  weights # weight matrix of assays / factors
  
  # for each view / assay, multiply the test assay data by the MOFA weight matrix
  test_data  = lapply(1:length(assays), function(x){
    d = t(test_data[[assays[[x]]]]) %*% MOFAobject@expectations$W[[assays[[x]]]]   # to produce a sample x factor matrix  
    d.w = lapply(1:ncol(d), function(y){
      d[,y] * weights[assays[[x]], y]
      
    })
    d.w = do.call(cbind, d.w)
  })
  
  # we now have factors for each sample for each view. we would prefer one factor for each sample, representing all views / assay 
  
  # sum weighted factors to give us a samples x factors matrix combining all assays
  test_data_sum = Reduce('+', test_data)
  dim(test_data_sum)
  test_data_sum[1:100, 1:15]
  
  test_data_sum = as.data.frame(test_data_sum)
  colnames(test_data_sum) = paste('Factor', 1:ncol(test_data_sum), sep='')
  return(test_data_sum)
}


