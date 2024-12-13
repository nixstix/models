source('../scripts/libs.R')
source('../scripts/run_glmnet.R')
detach("package:MOFA2", unload=TRUE) # need to detach MOFA in order to use the predict function from glmnet

load('../data/MOFA_models/predictors.RData')
load('../data/MOFA_models/test.RData')

# test predictions on 2022 data
# -----------------------------
# -----------------------------

# otehr to try
# features = c( 'Day0.Meta.age', 'Day0.TcellPol.pol', 'Day0.TcellActivation.PT')

task = 'IgG_PT'

## baseline model
## --------------

# use demographic data, as well as assay data that should be highly predictive of the prediction task


features = c( 'Day0.Meta.age', 'Day0.TcellPol.pol')
model1 = run_glmnet(p = features, task.name = task, alpha=1)
model1$model_coef
model1$model_cor

features = c('Day0.IgG_PT', 'Day0.Meta.age', 'Day0.TcellPol.pol')
model2 = run_glmnet(p = features, task.name = task, alpha=1)
model2$model_coef
model2$model_cor

features = c('Day0.IgG_PT', 'Day0.Meta.age', 'Day0.TcellPol.pol', 'Day0.TcellActivation.PT')
model3 = run_glmnet(p = features, task.name = task, alpha=1)
model3$model_coef
model3$model_cor

features = c('Day0.IgG_PT', 'Day0.TcellPol.pol', 'Day0.TcellActivation.PT')
model4 = run_glmnet(p = features, task.name = task, alpha=1)
model4$model_coef
model4$model_cor


features = c( 'Day0.TcellPol.pol', 'Day0.TcellActivation.PT')
model5 = run_glmnet(p = features, task.name = task, alpha=1)
model5$model_coef
model5$model_cor

features = c( 'Day0.TcellPol.pol', 'Day0.TcellActivation.PT', 'Day0.Meta.age')
model6 = run_glmnet(p = features, task.name = task, alpha=1)
model6$model_coef
model6$model_cor

save(model1, model2, model3, model4, model5,model6, file='../data/regression_models/regression_models_IgG_PT_noInt.RData')

model1$model_cor
model2$model_cor
model3$model_cor
model4$model_cor
model5$model_cor
model6$model_cor