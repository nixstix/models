The original pipeline is specified under `MOFA_7assays`. Any changes to this pipeline are documented below.

## MOFA/MOFA:

In the first chunk we have changed the root directory
In the first chunk we have changed the number of assays from 7 to 4

```
knitr::opts_knit$set(root.dir = '/mnt/BioAdHoc/Groups/Peters/nthrupp/CMI-PB3/MOFA_4assays/MOFA/')
nr_factors = 15

assays = c("geneExp",  "ab"  ,"cytokineL", "cell_freq")
```

## MOFA/prep_testData.R

At line 14 we have changed the number of assays from 7 to 4

```
assays = c('geneExp', 'ab',  'cytokineL', 'cell_freq')
```


