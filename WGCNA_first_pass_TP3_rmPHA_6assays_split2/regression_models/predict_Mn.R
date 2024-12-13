source('../scripts/libs.R')
source('../scripts/run_glmnet.R')
detach("package:MOFA2", unload=TRUE) # need to detach MOFA in order to use the predict function from glmnet

load('../data/MOFA_models/predictors.RData')
load('../data/MOFA_models/test.RData')

# test predictions on 2022 data
# -----------------------------
# -----------------------------

task = 'Freq.Monocytes'

## baseline model
## --------------

# use demographic data, as well as assay data that should be highly predictive of the prediction task


features = c("Day0.Freq.Monocytes", 'Day0.Meta.age')
model1 = run_glmnet(p = features, task.name = task, alpha=1)
model1$model_coef
model1$model_cor

## integration model
## -----------------

# use assay values from Day 0, as this has already been shown to be the most predictive,
# as well as factors generated from MOFA integration

# build model on 2021 data

f= grep('Day0.Factor', colnames(predictors.baseline), value=T)

# test all factors
p.int = c(f)
model2 = run_glmnet(p = p.int, task.name = task, alpha=1)
model2$model_cor
model2$model_coef

# test all factors plus baseline features
p.int = c(f, features)
model3 = run_glmnet(p = p.int, task.name = task, alpha=1)
model3$model_cor
model3$model_coef

# select top features from above

to_test = model2$model_coef
to_test = to_test[grep('Factor', names(to_test))]
to_test = to_test[order(abs(to_test), decreasing = T)]
to_test

models = lapply(1:length(to_test), function(x){
  to_test = to_test[1:x]
  p.int = c(features, names(to_test))
  model = run_glmnet(p = p.int, task.name = task, alpha=1)
  model$model_cor
  model$model_coef  
  return(model)
})

cors =  list()
for(i in 1:length(models)){
  cors[i] = models[[i]]$model_cor$estimate
}
cors = unlist(cors)
max(abs(cors))

models_alpha = lapply(1:length(models), function(i){
  if (abs(models[[i]]$model_cor$estimate) == max(abs(cors))){
    p.int = names(models[[i]]$model_coef)
    p.int = p.int[!p.int %in% c('(Intercept)')]
    model = lapply(c(0.0, 0.25, 0.5, 0.75, 1.0), function(x){
      run_glmnet(p = p.int, task.name = task, alpha=x)      
    })
  }  
})


save(model1, model2, model3, models, models_alpha, file='../data/regression_models/regression_models_Mn.RData')

if(exists('cors')){
  print(cors)
}

