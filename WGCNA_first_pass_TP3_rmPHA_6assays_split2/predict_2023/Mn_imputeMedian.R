# predict

source('../scripts/libs.R')

load('../data/predictions/challenge.RData')
load('../data/regression_models/regression_models_Mn.RData')

# which model:
mod = models[[4]]$model

predictors = rownames(coef(mod))
predictors = predictors[!predictors %in% '(Intercept)']

new_data = challenge_data.baseline[, c(predictors)]
dim(new_data)


new_data_i = data.frame(apply(new_data, 2, function(x){
  m = median(x, na.rm = T)
  x[is.na(x)] = m
  return(x)
})
)

new_data = na.omit(new_data_i)
dim(new_data_i)

preds = data.frame(predict(mod, newx=as.matrix(new_data_i), s='lambda.min'))

preds$rnk = rank(preds$lambda.min)
preds$subject_id = rownames(preds)
preds

all_subj = data.frame(subject_id=as.character(challenge_data.baseline$Day0.Meta.subject_id))

res = left_join(all_subj, preds)
res

write_tsv(x = res, file='../data/predictions/Mn_i.tsv', col_names = T)
