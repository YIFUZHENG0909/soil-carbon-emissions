# Install plspm package (uncomment if not installed)
# install.packages('plspm')
# devtools::install_github('gastonstat/plspm')

# Load plspm package
library(plspm)

# Set working directory
setwd("D:/Carbon Flux Data/data_upload/Fig5")

# Read the data
dat <- read.delim('data.txt', sep = '\t')

# Specify latent variables and their relationships in R as a list
dat_blocks <- list(
  Factory = 'Distance',
  soil = c('TC', 'CN'),
  Actinobacteria = 'Actinobacteria',
  shannon = 'Shannon',
  Gene_CD = c('K01179', 'K01178', 'K01176', 'K01209', 'K01191', 'K06113', 'K01205'),
  cf = "cf"
)
dat_blocks

# Describe relationships between latent variables using a 0-1 matrix,
# where 0 means no relationship and 1 means there is a relationship
Factory <- c(0,0,0,0,0,0)
soil <- c(1,0,0,0,0,0)
Actinobacteria <- c(1,1,0,0,0,0)
shannon <- c(0,1,1,0,0,0)
Gene_CD <- c(0,1,1,1,0,0)
cf <- c(0,0,0,0,1,0)

# Combine the paths
dat_path <- rbind(Factory, soil, Actinobacteria, shannon, Gene_CD, cf)
colnames(dat_path) <- rownames(dat_path)
dat_path

# Specify causal relationships, where A represents the columns as predictors and B represents the rows as outcomes
dat_modes <- rep('A', 6)
dat_modes

# A simple PLS-PM model; for more details, refer to ?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)

# Review the path coefficients and associated statistical information
dat_pls$path_coefs
dat_pls$inner_model

# Plot the causal relationship diagram (inner plot)
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

# View the status of endogenous and exogenous latent variables
dat_pls$inner_summary

# View the effect between variables
dat_pls$effects

# View the relationship between observed variables and latent variables; can also be visualized using outerplot()
dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

# Goodness-of-fit value to evaluate the model's fit
dat_pls$gof

# View latent variable scores, which are the standardized values of the latent variables
dat_pls$scores

