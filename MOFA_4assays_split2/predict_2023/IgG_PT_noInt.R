# predict

source('../scripts/libs.R')

load('../data/predictions/challenge.RData')
load('../data/regression_models/regression_models_IgG_PT_noInt.RData')

# which model:
model1$model_cor
mod = model1$model

predictors = rownames(coef(mod))
predictors = predictors[!predictors %in% '(Intercept)']

new_data = challenge_data.baseline[, c(predictors)]
dim(new_data)

new_data = na.omit(new_data)
dim(new_data)

preds = data.frame(predict(mod, newx=as.matrix(new_data), s='lambda.min'))

preds$rnk = rank(preds$lambda.min)
preds$subject_id = rownames(preds)
preds

all_subj = data.frame(subject_id=as.character(challenge_data.baseline$Day0.Meta.subject_id))

res = left_join(all_subj, preds)
res

write_tsv(x = res, file='../data/predictions/IgG_PT_noInt.tsv', col_names = T)
