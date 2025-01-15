library(randomForest)
##read table
otu <- read.table("key_gene.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
cf <- read.delim('carbon_flux(normalize to 0-1).txt', row.names = 1)
otu <- data.frame(t(otu))
otu <- otu[rownames(cf), ]
otu <- cbind(otu, cf)


# Divide the total data set into training sets (70%) and test sets (30%)
set.seed(123)
train <- sample(nrow(otu), nrow(otu)*0.7)
otu_train <- otu[train, ]
otu_test <- otu[-train, ]

# Random forest computation
set.seed(123)
otu_train.forest <- randomForest(cf~., data = otu_train, importance = TRUE, ntree = 1000)
#############################################################################################################
# Use the training set to see the prediction accuracy
plant_predict <- predict(otu_train.forest, otu_train)

plot(otu_train$cf, plant_predict, main = 'training',
     xlab = 'cf', ylab = 'Predict')
abline(1, 1)

# Use test sets to evaluate predictive performance
plant_predict <- predict(otu_train.forest, otu_test)

plot(otu_test$cf, plant_predict, main = 'test',
     xlab = 'cf', ylab = 'Predict')
abline(1, 1)

# Extract the importance of each variable for sample classification
imp = data.frame(importance(otu_train.forest),MDA.p = otu_train.forest$importanceSD[4])
head(imp) 

# Arrange variables in descending order of MeanDecreaseGini importance
library(dplyr)
imp = arrange(imp,desc(X.IncMSE)) 

## Output importance ranking results
write.csv(imp,"importance.csv",quote = FALSE)

## Cross-validation AIDS in evaluating the selection of a specific number of OTUs
# repeat ten fold cross verification
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu_train[-ncol(otu_train)], otu_train$cf, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv

# Extract the verification result plot
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 30)

library(ggplot2)

plot <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of Gunes', y = 'Cross-validation error')
plot
ggsave("cross.pdf", plot, width = 3, height = 3)