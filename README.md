# MicrobiomeFRM
MicrobiomeFRM is an R library for analyzing the Beta-diversity in microbiome data, using a semiparametric functional response model (FRM)

## Installation
```{r}
devtools::install_github('y1zhong/MicrobiomeFRM')
```
## Usage

```R
beta.dev(raw.OTU.datframe, 'bray') # returns 'Bray-Curtis" type of Beta-diversity with raw OTU counts 'raw.OTU.datframe'
pairwise.dist(df.cts, "euclidean") # returns 'Euclidean" type of between-subject difference with continuous type of variable in 'df.cts'
ugeecov_cate_cont_exp2_new(d.mat, df.cat, dist.list) # returns parameter estimates and asymptotic variances to analyze the association between the Beta-diversity 'd.mat' and categorical explanatory variables 'df.cat', as well as a list of between-subject attributes for continuous explanatory variables 'dist.list'
```
