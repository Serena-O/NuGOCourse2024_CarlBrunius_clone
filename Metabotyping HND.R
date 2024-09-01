# Check that you have the right R packages installed - and install them if not!
if (!require("remotes", quietly = TRUE)) install.packages('remotes')
if (!require("vioplot", quietly = TRUE)) install.packages('vioplot')
if (!require("factoextra", quietly = TRUE)) install.packages('factoextra')
if (!require("MUVR2", quietly = TRUE)) remotes::install_github('MetaboComp/MUVR2')
if (!require("rdCV", quietly = TRUE)) remotes::install_gitlab('CarlBrunius/rdCV')

# Clean the slate
rm(list = ls())

# Load the "Healthy Nordic Diet" synthetic data set from the triplot package
load(file = 'HealthyNordicDiet.rda')

# Extract metabolite data
MB <- HealthyNordicDiet$MetaboliteData
View(MB)

# Make a PCA of the metabolite data
pca <- prcomp(MB, center = TRUE, scale = TRUE)

# PCA biplot
rdCV::biplotPCA(pca, comps = 1:2, labPlSc = FALSE)

# Extract baseline data
BL <- HealthyNordicDiet$BaselineData
View(BL)
BL[c(1, 501),]
BL[c(5, 505),]

# PCA biplot colored by baseline characteristics - continuous
rdCV::biplotPCA(pca, comps = 1:2, xCol = BL$Energy, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, xCol = BL$Age, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, xCol = BL$BMI, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, xCol = BL$FastingGlucose, labPlSc = FALSE)

# PCA biplot colored by baseline characteristics - factors
rdCV::biplotPCA(pca, comps = 1:2, colSc = BL$Group, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, colSc = BL$Gender, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, colSc = BL$Smoking, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, colSc = BL$PhysicalActivityIndex, labPlSc = FALSE)
rdCV::biplotPCA(pca, comps = 1:2, colSc = BL$Education, labPlSc = FALSE)





# Choose variables to do metabotyping (clustering)
# Cluster based on age, BMI and fasting glucose
clusterData <- scale(BL[, c(6, 7, 11)])
clusterData <- as.data.frame(clusterData)
head(clusterData)

# K-means clustering with silhouette statistic
library(factoextra)
fviz_nbclust(clusterData, kmeans, method='silhouette')
# Silhouette scores indicate 2 clusters
k <- kmeans(clusterData, 2)
metabotype <- as.factor(k$cluster)
table(metabotype)

# Examine how metabotypes relate to the clustering variables
library(vioplot)
par(mfrow = c(3, 1), mar = c(4, 4, 0, 0) + 0.5)
vioplot(clusterData$Age ~ metabotype, col = c('green', 'red'), ylab = 'Age')
(ttest <- t.test(clusterData$Age[metabotype == 1], clusterData$Age[metabotype == 2]))
legend('topleft', paste('p =', sprintf("%.3f", ttest$p.value)), bty = 'n')
vioplot(clusterData$BMI ~ metabotype, col = c('green', 'red'), ylab = 'BMI')
(ttest <- t.test(clusterData$BMI[metabotype == 1], clusterData$BMI[metabotype == 2]))
legend('topleft', paste('p =', sprintf("%.3f", ttest$p.value)), bty = 'n')
vioplot(clusterData$FastingGlucose ~ metabotype, col = c('green', 'red'), ylab = 'Fasting glucose')
(ttest <- t.test(clusterData$FastingGlucose[metabotype == 1], clusterData$FastingGlucose[metabotype == 2]))
legend('topleft', paste('p =', sprintf("%.3f", ttest$p.value)), bty = 'n')


# PCA biplot colored by metabotype
par(mfrow = c(1, 1), mar = c(4, 4, 0, 0) + 0.5)
rdCV::biplotPCA(pca, comps = 1:2, colSc = metabotype, labPlSc = FALSE)






# Some quick analyses to investigate metabotype attributes
par(mfrow = c(1, 1), mar = c(4, 4, 0, 0) + 0.5)
vioplot(BL$Energy ~ metabotype, col = c('green', 'red'), ylab = 'Total Energy') # A little higher energy intake in metabotype 1
ttest <- t.test(BL$Energy[metabotype == 1], BL$Energy[metabotype == 2]) 
legend('topleft', paste('p =', sprintf("%.3f", ttest$p.value)), bty = 'n')

makeChisq <- function(variable) {
  variableName <- deparse(substitute(variable))
  cat('Chi2 test for variable:', variableName, '\n\n')
  c <- chisq.test(variable, metabotype)
  names(dimnames(c$observed))[1] <- substring(variableName, 4)
  print(c$observed)
  cat('\np-value\n')
  cat(c$p.value)
}

makeChisq(BL$Group) # T2D strongly associated to metabotype
makeChisq(BL$Gender) # Gender not associated to metabotype
makeChisq(BL$Smoking) # Smoking not associated to metabotype
makeChisq(BL$PhysicalActivityIndex) # PA borderline associated to metabotype
makeChisq(BL$Education) # Education associated to metabotype





# Associate metabotype to metabolite profile
# Call in packages
library(MUVR2)
library(doParallel)
# Set up for parallel processing
nCore <- detectCores() - 1
cl <- makeCluster(nCore)
registerDoParallel(cl)
# Actual model
HND_model <- MUVR2(X = MB, 
                   Y = metabotype, 
                   method = 'RF', 
                   nRep = 3 * nCore, 
                   varRatio = 0.75)
# Stop parallel processing
stopCluster(cl)

# Examine model
plotVAL(HND_model)
HND_model$fitMetric
HND_model$nVar
confusionMatrix(HND_model, model = 'min')
head(cbind(metabotype, HND_model$yClass), 15)

# Examine variables of interest
plotVIRank(HND_model, model = 'min', add_blank = 14)
getVIRank(HND_model, model = 'min')







# Extract Variables-of-Interest
VoI <- subset(MB, select = getVIRank(HND_model, model = 'min')$name)
head(VoI)

# Investigate VoIs in relation to metabotype
par(mfrow = c(5, 4), mar = c(4, 4, 0, 0) + 0.5)
for(metab in 1:ncol(VoI)) {
  vioplot(VoI[,metab][[1]] ~ metabotype, col = c('green', 'red'), ylab = colnames(VoI)[metab])
  ttest <- t.test(VoI[metabotype == 1, metab], VoI[metabotype == 2, metab])
  legend('topleft', paste('p =', sprintf("%.3f", ttest$p.value)), bty = 'n')
}






# Associate metabotypes to diet
diet <- HealthyNordicDiet$FoodData
View(diet)
par(mfrow = c(5, 4), 
    mar = c(4, 4, 0, 0) + 0.5)
for(dietItem in 1:ncol(diet)) {
  vioplot(diet[,dietItem][[1]] ~ metabotype, col = c('green', 'red'), ylab = colnames(diet)[dietItem])
  wtest <- wilcox.test(diet[metabotype == 1, dietItem][[1]], diet[metabotype == 2, dietItem][[1]])
  legend('topleft', paste('p =', sprintf("%.3f", wtest$p.value)), bty = 'n')
}
