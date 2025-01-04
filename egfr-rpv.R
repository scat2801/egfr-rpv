#sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-devlibcairo2-dev libgif-dev
#sudo apt-get install libgfortran4
#sudo apt-get install libcairo2-dev libjpeg-dev libgif-dev

install.packages('sysfonts')
install.packages('showtext')

install.packages(c("FactoMineR", "alookr", "Boruta", "factoextra", "ggofortify", "randomForest", "praznik", "mboost", "mltools", "VIM", "pheatmap", "neuroCombat", "praznik", "caretEnsemble"))

install.packages("randomForest")

devtools::install_github("RomeroBarata/bimba")

install.packages("gdtools", type = "source")
devtools::install_github("dreamRs/esquisse")

install.packages("alookr")
install.packages("ggfortify")
install.packages("devtools", force = TRUE)

library(devtools)
install_github("jfortin1/neuroCombat_Rpackage", force = TRUE)

install.packages("igraph")
install.packages("pbkrtest")

install.packages('BiocManager')
BiocManager::install('RegParallel', force = TRUE)

install.packages("alookr") 
install.packages("gtable")
install.packages("caTools")
install.packages("glmnet")
install.packages("smotefamily")
install.packages("gtsummary")
install.packages("cards")
install.packages("cardx")
install.packages("fpp3")

library(cards)
library(cardx)
library(fpp3)
library(gtsummary)
library(smotefamily)
library(caTools)
library (gtable)
library("FactoMineR")
library("factoextra")
library(e1071)
library(survival)
library(RegParallel)
library(bimba)
library(purrr)
library(glmnet)
library(dplyr)  
library(magrittr)
library(knitr)
library(survMisc)
library(ggplot2)
library(ggfortify)
library(data.table)
library(caret)
library(boot)
library(class)
library(randomForest)
library(ROCR)
library(pROC)
library(Boruta)
library(praznik)
library(xgboost)
library(mboost)
library(RColorBrewer)
library(mltools)
library(VIM)
library(readxl)
library(mice)
library(corrplot)
library(pheatmap)
library(alookr)
library(devtools)
library(Boruta)
library(neuroCombat)
library(praznik)
library(caretEnsemble)
library(glmnet)
library(magrittr) # needs to be run every time you start R and want to use %>%

'%ni%' <- Negate('%in%')

#######################Import Main Data#############

datapath <- "/home/mitchchen/git/LCHisto/radiomics_opsheet.xlsx"

lcdata <- read_excel(datapath, sheet = "radiomics-i")
lcdata$Performance <- as.numeric (lcdata$Performance)

#Rename AUC columns
names(lcdata)[names(lcdata) == "AUC-CSH"] <- "AUC_CSH"
names(lcdata)[names(lcdata) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(lcdata)[names(lcdata) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(lcdata)[names(lcdata) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(lcdata)[names(lcdata) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(lcdata)[names(lcdata) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(lcdata)[names(lcdata) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(lcdata)[names(lcdata) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(lcdata)[names(lcdata) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(lcdata)[names(lcdata) == "N voxels"] <- "N_voxels"

#Start of the feature columns
st_col <- 70

for (i in 1:nrow(lcdata)){
  lcdata$Stage[i] <- substring(lcdata$Stage[i], 1, 1)
  lcdata$N[i] <- substring(lcdata$N[i], 1, 1)
  lcdata$M[i] <- substring(lcdata$M[i], 1, 1)
}

#Import shell data
lcdata_shell <- read_excel(datapath, sheet = "radiomics-s")

names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH"] <- "AUC_CSH"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(lcdata_shell)[names(lcdata_shell) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(lcdata_shell)[names(lcdata_shell) == "N voxels"] <- "N_voxels"

#Import patch data
lcdata_patch <- read_excel(datapath, 
                           sheet = "radiomics-p")

names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH"] <- "AUC_CSH"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(lcdata_patch)[names(lcdata_patch) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(lcdata_patch)[names(lcdata_patch) == "N voxels"] <- "N_voxels"

##############Combine the three masks#################
lcdata_interim <- merge(lcdata, lcdata_shell, by = "Patient.ID", all=TRUE)
lcdata_allfeatures <- merge(lcdata_interim, lcdata_patch, by = "Patient.ID", all=TRUE)

rm(lcdata_interim, lcdata_shell, lcdata_patch)

#x: lesion, .y: shell, .: patch

#Use for unsupervised clustering
lcdata <- subset(lcdata_allfeatures, Contrast == 1)
lcdata <- subset(lcdata, Thickness == 0)

#Divide into surgical vs non-surgical
surgical <- subset(lcdata_allfeatures, Csurg == 1)

#Space and time separation and external validation
SMH <- subset(lcdata_allfeatures, Space == 1)
HH <- subset(lcdata_allfeatures, Space == 2)
CXH <- subset(lcdata_allfeatures, Space == 3)
OH <- subset(lcdata_allfeatures, Space == 4)
Pre<- subset(lcdata_allfeatures, Time == 0)
Post<- subset(lcdata_allfeatures, Time == 1)

####################Feature Class import########################################
#Covariate import
covariates <- c(colnames(lcdata))
covariates <- covariates[st_col:ncol(lcdata)]

###################Import external validation sets############

#Main lesion data

NRdata <- read_excel(datapath, 
                     sheet = "TV-i")
names(NRdata)[names(NRdata) == "AUC-CSH"] <- "AUC_CSH"
names(NRdata)[names(NRdata) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(NRdata)[names(NRdata) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(NRdata)[names(NRdata) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(NRdata)[names(NRdata) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(NRdata)[names(NRdata) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(NRdata)[names(NRdata) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(NRdata)[names(NRdata) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(NRdata)[names(NRdata) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(NRdata)[names(NRdata) == "N voxels"] <- "N_voxels"

NRdata_shell <- read_excel(datapath, 
                     sheet = "TV-s")
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH"] <- "AUC_CSH"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(NRdata_shell)[names(NRdata_shell) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(NRdata_shell)[names(NRdata_shell) == "N voxels"] <- "N_voxels"

NRdata_patch <- read_excel(datapath, 
                           sheet = "TV-p")
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH"] <- "AUC_CSH"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_LHL"] <- "AUC_CSH_LHL"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_LHH"] <- "AUC_CSH_LHH"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_LLH"] <- "AUC_CSH_LLH"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_LLL"] <- "AUC_CSH_LLL"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_HLL"] <- "AUC_CSH_HLL"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_HLH"] <- "AUC_CSH_HLH"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_HHL"] <- "AUC_CSH_HHL"
names(NRdata_patch)[names(NRdata_patch) == "AUC-CSH_HHH"] <- "AUC_CSH_HHH"
names(NRdata_patch)[names(NRdata_patch) == "N voxels"] <- "N_voxels"

NRdata_interim <- merge(NRdata, NRdata_shell, by = "Patient.ID", all=TRUE)
NRdata_allfeatures <- merge(NRdata_interim, NRdata_patch, by = "Patient.ID", all=TRUE)

rm(NRdata, NRdata_interim, NRdata_shell, NRdata_patch)

data <- subset(lcdata_allfeatures, !is.na(EGFR))
rm(lcdata_allfeatures)

extdata <- subset(NRdata_allfeatures, !is.na(EGFR))

st_col_ext <- 25

#EGFR <- data$EGFR
rm(NRdata_allfeatures)

data$Age[65] <- 66

#Summary statistics
mean (as.numeric(extdata$Age))
sd (as.numeric(extdata$Age))
max (as.numeric(extdata$Age))
min (as.numeric(extdata$Age))

tbl <- table(extdata$Gender)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$Stage)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$N)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$M)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$Histology)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$PDL1)
cbind(tbl,prop.table(tbl))

tbl <- table(extdata$EGFR)
cbind(tbl,prop.table(tbl))

wilcox.test(x = data$Age, y = as.numeric(extdata$Age),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


wilcox.test(x = as.numeric(data$EGFR), y = as.numeric(extdata$EGFR),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


#############
### SPLIT DATA INTO TRAINING AND EXTERNAL VALIDATION 
set.seed(43)

data$His <- data$EGFRMut

#External Validation Data

for (i in 1:nrow(extdata)) {
  if (extdata$EGFR[i] == 1) {
    if((extdata$KRAS[i] == 0)&(extdata$ALK[i] == 0)){
      extdata$EGFRMut[i] <- 1
    }
  }
}

split = sample.split(data$His, SplitRatio = 0.8)
train_set_full <- subset(data, split == T)
val_set_full <- subset(data, split == F)
train_set_full <- as.data.frame(train_set_full)
val_set_full <- as.data.frame(val_set_full)

#st_col <- 70
head(train_set_full[st_col:ncol(data)])

#move last columns first
train_set_full <- train_set_full %>%
  dplyr::select(EGFRMut, everything())

train_set_full <- train_set_full %>%
  dplyr::select(His, everything())

val_set_full <- val_set_full %>%
  dplyr::select(EGFRMut, everything())

val_set_full <- val_set_full %>%
  dplyr::select(His, everything())

extdata <- extdata %>%
  dplyr::select(EGFRMut, everything())

st_col <- 72

#train_set_full <- train_set_full_original 

#train_set_full <- train_set_full[st_col:ncol(data)]
#SCALE DATA

r_train_mean <- colMeans(train_set_full[st_col:ncol(train_set_full)]) ## calulate the mean of the numeric features, note not the categorical features
r_train_std <- sapply(train_set_full[st_col:ncol(train_set_full)], sd, na.rm = TRUE) ## Standard deviation
r_scaled_train = as.data.frame(scale(train_set_full[st_col:ncol(train_set_full)], center=r_train_mean, scale= r_train_std))

r_val_mean <- colMeans(val_set_full[st_col:ncol(val_set_full)]) 
r_val_std <- sapply(val_set_full[st_col:ncol(val_set_full)], sd, na.rm = TRUE) 
r_scaled_val  = as.data.frame(scale(val_set_full[st_col:ncol(val_set_full)], center= r_val_mean, scale= r_val_std)) # scale the external validation using training data set statistics

st_col_ext <- 26
r_ext_mean <- colMeans(extdata[st_col_ext:ncol(extdata)]) 
r_ext_std <- sapply(extdata[st_col_ext:ncol(extdata)], sd, na.rm = TRUE) 
r_scaled_ext  = as.data.frame(scale(extdata[st_col_ext:ncol(extdata)], center= r_ext_mean, scale= r_ext_std)) 

### Feature reduction/dimensionality reduction steps are either:
# unsupervised - (1) PCA on the full radiomics feature set
# supervised - removing highly correlated features, FDR step and then: (2) LASSO/ (3) E-NET/ (4) RFE/ (5) BORUTA

x <- as.matrix(r_scaled_train)  ### Here we are selecting just the input variables, not the outcome variable, as PCA is an unsupervised technique

x <- x[ , which(apply(x, 2, var) != 0)]
r_scaled_train <- r_scaled_train[ , which(apply(r_scaled_train, 2, var) != 0)]

### Here we undertake principal component analysis for dimensionality reduction
prin_comp <- prcomp(x, scale. = F)  
summary (prin_comp)

fviz_pca_var(prin_comp, col.var = "black")
fviz_cos2(prin_comp, choice = "var", axes = 1:2)

# calculate the std deviation of the PCs and then compute the variance
std_dev <- prin_comp$sdev
pr_var <- std_dev^2

# compute the proportion of the variance explained by each PC
prop_varex <- pr_var/sum(pr_var)

# create a scree plot to look at how many PCs to select - "zoom" in to the first x number of PCs rather than looking at all
zoom <- prop_varex[1:50]
plot(zoom, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

# A cumulative scree plot also helps to visualize - this shows that the first ~X PCs explain approx 90% of the variance
plot(cumsum(zoom), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

res.pca <- PCA(x, graph = T, scale.unit = T, ncp = 5)
#print(res.pca)
eig.val <- get_eigenvalue(res.pca)

fviz_eig(res.pca, addlabels = T, ylim = c(0,50))

vari <- get_pca_var(res.pca)

# Contributions of variables to PC1&2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 15)

### REMEMBER TO CHANGE THE NUMBER OF PCs SELECTED HERE - ACCOUNT FOR 80% of VARIANCE OR 15 PCs - whichever is lowest

### Now to make new feature sets for modeling which are the PCs
train.pca<- as.data.frame(prin_comp$x)
train.pca<- train.pca[,1:15] #### HERE YOU SELECT THE NUMBER OF PCs THAT ACCOUNT FOR e.g. 90% of the variance
train.pca$His <- train_set_full$His

val.pca <- as.data.frame(predict(prin_comp, newdata = r_scaled_val))
val.pca <- val.pca[,1:15] #### AS ABOVE - HERE YOU SELECT THE NUMBER OF PCs THAT ACCOUNT FOR e.g. 90% of the variance
val.pca$His <- val_set_full$His

ev.pca <- as.data.frame(predict(prin_comp, newdata = r_scaled_ext))
ev.pca <- ev.pca[,1:15] #### AS ABOVE - HERE YOU SELECT THE NUMBER OF PCs THAT ACCOUNT FOR e.g. 90% of the variance
ev.pca$His <- extdata$His

autoplot(prin_comp)

autoplot(prin_comp, data = train.pca, colour = 'His') + theme(title = element_text(size = 15), legend.text = element_text(size = 12)) + theme_bw()

### NOW FOR THE SUPERVISED APPROACHES - FOR THIS WE FIRST REMOVE HIGHLY CORRELATED FEATURES AND APPLY THE FDR STEP

### Now to look at just the radiomic features and remove any that are highly correlated
tmp <- cor(r_scaled_train)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
# # Above two commands can be replaced with
#  tmp[!lower.tri(tmp)] <- 0
# #
harm_rad_train2 <- as.data.frame(as.matrix(r_scaled_train)[,!apply(tmp,2,function(x) any(x > 0.90))])

head(harm_rad_train2)
ncol(harm_rad_train2)

harm_rad_train2$His <- as.factor(train_set_full$His)

harm_rad_train2 <- harm_rad_train2 %>%
  dplyr::select(His, everything())

#harm_rad_train2 <- harm_rad_train2[-1]

### Now we undertake the FDR step by taking the univariate LR for each variable in the training set. We extract the list of p-values and apply BH correction to account for multiple testing. Then we create a list of significant values (number of values differs according to how strict a p value threshold is used). This can take several minutes.
#res = tbl_uvregression(data = harm_rad_train2, method = glm, method.args = list(family = binomial), y = His)

#res_sig = res$meta_data
#res_sig$p = res[["table_body"]][["p.value"]]
#res_sig = res_sig[,-c(2,3,4)]
#res_sig$p_adjusted = p.adjust(res_sig$p, method = "BH", n = length(res_sig$p))

#final = subset(res_sig, p < 0.05)
#uni_LR_variables <- final$variable

#length(uni_LR_variables)

#################Alternative FDR step######################
covariates <- c(colnames(harm_rad_train2))
covariates <- covariates[2:ncol(harm_rad_train2)]

#univ_formulas <- sapply(covariates, function(x) as.formula(paste('`His`~', x)))
#univ_models <- lapply(univ_formulas, function(x){glm(x, data = harm_rad_train2, family = binomial(link = logit))})

#univ_results <- lapply(univ_models,
#                       function(x){ 
#                         x <- summary(x)
#                         p.value<-signif(x$coef[,4], digits=5)
#                         beta<-signif(x$coef[1], digits=5);#coefficient beta
#                         esta <-signif(x$coef[2], digits=5);#exp(beta)
#                         res <- c(beta, esta, p.value)
#                         names(res)<-c("beta", "estimate", "p.value")
#                         return(res)
#                       })
#res <- t(as.data.frame(univ_results, check.names = FALSE))
############

#Preferred FDR Method
set.seed(42)
univ_formulas <- sapply(covariates,
               function(x) {
                 tmp <- try(coef(summary(glm(as.formula(paste("His", x, sep = "~")),
                                             family=binomial, data = harm_rad_train2)))[2, ], TRUE)
                 if (class(tmp) == "try-error") NULL else tmp})

# 
#rm(res)
res <-  as.data.frame(t(univ_formulas))

#significance 0.10
final = subset(res, res[,4] < 0.05)
uni_LR_variables<-rownames(final)
length(uni_LR_variables)

### Now we subset the training and validation sets to use the significant variables
harm_rad_train2 <- harm_rad_train2[, uni_LR_variables]
harm_rad_val2 <- as.data.frame(r_scaled_val[,uni_LR_variables])

extHis <- extdata$EGFRMut
extdata <- as.data.frame(r_scaled_ext[,colnames(harm_rad_train2)])
extdata$His = extHis


x <- as.matrix(harm_rad_train2)
y = as.numeric(train_set_full$His)

### Now we we will perform a LASSO (10-fold cross validation on the training data)
set.seed(123)

na_index <- is.na(y)

x <- x[!na_index, ]
y <- y[!na_index]

fit <- glmnet(x, y, family = "binomial", alpha = 1, standardize=TRUE)

cv <- cv.glmnet(x, y, family="binomial",alpha=1, standardize=TRUE)
fit <- glmnet(x, y, family="binomial",alpha=1,lambda=cv$lambda.min,  standardize=TRUE)

#dev.copy(jpeg, 'lambda_lasso_Plot.jpg')
plot(cv)
#dev.off()
pred <- predict(fit, x)

c<-coef(fit,s='lambda.min') ### can change this to  lambda.1se
inds<-which(c!=0)
variables<-row.names(c)[inds]

'%ni%' <- Negate('%in%')
Lasso.vars <-variables[variables %ni% '(Intercept)']
print(length(Lasso.vars))
Lasso.vars ## This will print the list of variables here... 

### Now we perform ELASTIC-NET (10-fold cross validation on the training data)
set.seed(123)
fit <- glmnet(x, y, family="binomial",alpha=0.5, standardize=TRUE) ## can always loop throguh alpha but 0.5 is okay for a stating point 
#dev.copy(jpeg, 'Lasso_Trace_Plot.jpg')
plot(fit, label=)
#dev.off()
cv <- cv.glmnet(x, y, family="binomial",alpha=0.5,  standardize=TRUE)
fit <- glmnet(x, y, family="binomial",alpha=0.5,lambda=cv$lambda.min,  standardize=TRUE) ## can change to lambda.1se but min is best
#dev.copy(jpeg, 'lambda_lasso_Plot.jpg')
plot(cv)
#dev.off()
pred <- predict(fit, x)
c<-coef(fit,s='lambda.min') ### again can change this to  lambda.1se
inds<-which(c!=0)
variables<-row.names(c)[inds]
E.Net.vars <-variables[variables %ni% '(Intercept)']
print(length(E.Net.vars))
E.Net.vars ## print the variables here...

### Recursive Feature Elimination
x <- as.data.frame(x)
y <- as.factor(y)

set.seed(321)
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # repeated cv
                      repeats = 5, # number of repeats
                      number = 10) # number of folds

result_rfe1 <- rfe(x, y, sizes = c(1:116), rfeControl = control)

# Print the results
result_rfe1

# Print the selected features
predictors(result_rfe1)
rfe.vars <- predictors(result_rfe1)

# Print the results visually
ggplot(data = result_rfe1, metric = "Accuracy") + theme_bw()
ggplot(data = result_rfe1, metric = "Kappa") + theme_bw()


### Here I plot the top 10 features in order of importance as per RFE (can change to more than 10 if needed)
varimp_data <- data.frame(feature = row.names(varImp(result_rfe1))[1:20],
                          importance = varImp(result_rfe1)[1:20, 1])



rfe_importance <- ggplot(data = varimp_data,
                         aes(x = reorder(feature, importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance", title="Recursive Feature Elimination Results") +
  geom_text(aes(label = round(importance, 2)), hjust = 1.2, color="black", size=3) +
  theme_bw() + theme(legend.position = "none")

rfe_importance + coord_flip()

#5) Boruta
set.seed(123)
y <- as.factor(y)
boruta.train <- Boruta(x, y, doTrace = 2)

final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
bor.vars <- getSelectedAttributes(final.boruta, withTentative = F)

#6) Mutual Information

set.seed(123)

MIM.train <- MIM(x, y, k=10)
mim.vars <- MIM.train$selection
mim.vars

#7) Pearson
### Feature Selection with Pearson Correlation
outcome <- as.numeric(y)

pearson_var <- cor(x,outcome, method = "pearson")
pearson_var <- as.data.frame(abs(pearson_var))
pearson_var <- top_frac(pearson_var, 0.2)
p.vars <- row.names(pearson_var)
print(length(p.vars))
p.vars

#8) Spearman

### Feature Selection with Spearman Correlation
spearman_var <- cor(x,outcome, method = "spearman")
spearman_var <- as.data.frame(abs(spearman_var))
spearman_var <- top_frac(spearman_var, 0.2)
s.vars <- row.names(spearman_var)
print(length(s.vars))
s.vars

#9) Kendall
kendall_var <- cor(x,outcome, method = "kendall")
kendall_var <- as.data.frame(abs(kendall_var))
kendall_var <- top_frac(pearson_var, 0.2)
k.vars <- row.names(pearson_var)
print(length(k.vars))
k.vars

#Regression

pca = str("pca")
features_select <- list(p.vars, s.vars, k.vars, rownames(as.data.frame(mim.vars)), Lasso.vars, E.Net.vars, bor.vars, rfe.vars, pca)

#without pca
#features_select <- list(p.vars, s.vars, k.vars, rownames(as.data.frame(mim.vars)), Lasso.vars, E.Net.vars, bor.vars, rfe.vars)

res <- NULL
res_ext <- NULL

j <- 0

for (i in 1:length(features_select)){
  j = j + 1
  print(i)

#With PCA
  if (j == length(features_select)){ 
    ### Run this chunk in order to convert train2 and test2 into the PCA training and test data - to use PCA for the following models
    rad_train3 <- NULL
    rad_train3 = train.pca
    rad_val3 <- NULL
    rad_val3 = val.pca
    rad_ext <- NULL
    rad_ext <- ev.pca
    j <- 0
  }
  
  else if(j != length(features_select)){
    rad_train3 <- NULL
    rad_train3 <- harm_rad_train2[ ,unlist(features_select[i])]
    rad_val3 <- NULL
    rad_val3 <- harm_rad_val2[ ,unlist(features_select[i])]
    rad_ext <- NULL
    rad_ext <- extdata[ ,unlist(features_select[i])]
  }
  #rad_train3<- NULL
  #rad_train3 <- harm_rad_train2[ ,i] ### Use "Lasso.vars" or "E.Net.Vars" or "rfe.vars" or "bor.vars" etc to switch from one to the other
  rad_train3$His <- y
  
  #rad_val3<- NULL
  #rad_val3 <- harm_rad_val2[ ,i] ### Use "Lasso.vars" or "E.Net.Vars" or "rfe.vars" or "bor.vars" etc to switch from one to the other
  rad_val3$His <- val_set_full$His 
  
  synth_train <- rad_train3
  
  rad_ext$His <- extdata$His
  
  #synth_train <- ADASYN(rad_train3, perc_min = 50, k=5)
  #synth_val <- ADASYN(rad_val3, perc_min = 50, k=5)
  
  #rad_train3 <- synth_train
  ###rad_val3 <- synth_val

  ### 1) GENERALISED LINEAR MODEL
  ###Fitting logistic regression to the training set
   rad_train3$His <- as.factor(rad_train3$His)
   rad_val3$His <- as.factor(rad_val3$His)
   rad_ext$His <- as.factor(rad_ext$His)
   
   trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
   set.seed(42)
   classifier_glm <- train(form = His ~ .,
                           data = rad_train3,
                           method = "glm",
                           family=binomial,
                           trControl = trctrl,
                           tuneLength = 20
  )
   
  # #Predicting the validation set results
   glm_pred = predict(classifier_glm, newdata = rad_val3[-ncol(rad_val3)])
   glm_pred = factor(glm_pred, levels = c("0","1"), labels = c("0","1"))
   cm_glm = confusionMatrix(glm_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
   glm_pred = predict(classifier_glm, newdata = rad_ext[-ncol(rad_ext)])
   glm_pred = factor(glm_pred, levels = c("0","1"), labels = c("0","1"))
   cm_glm_ext = confusionMatrix(glm_pred, reference = as.factor(rad_ext$His), positive = "1") 
   
   
  ## 2 KNN
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  set.seed(42)
  
  ###Fitting knn to the training set
  classifier_knn2 = train(form = His ~ .,
                          data = rad_train3,
                          method = "knn",
                          trControl=trctrl,
                          tuneLength = 20)
  
  classifier_knn2$bestTune #Predicting the validation set results
  knn2_pred = predict(classifier_knn2, newdata = rad_val3[-ncol(rad_val3)])
  knn2_pred = factor(knn2_pred, levels = c("0","1"), labels = c("0","1"))
  cm_knn2 = confusionMatrix(knn2_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation
  
  knn2_pred = predict(classifier_knn2, newdata = rad_ext[-ncol(rad_ext)])
  knn2_pred = factor(knn2_pred, levels = c("0","1"), labels = c("0","1"))
  cm_knn2_ext = confusionMatrix(knn2_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation
  
  ### 3) LINEAR SUPPORT VECTOR MACHINE - with GridSearch using CARET (optimised hyperparameters for better results)
  set.seed(42)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  classifier_svm2 = train(form = His ~ .,
                          data = rad_train3,
                          method = "svmLinear",
                          trControl=trctrl,
                          tuneLength = 20
  )
  svm_pred = predict(classifier_svm2, newdata = rad_val3)
  svm_pred = factor(svm_pred, levels = c("0","1"), labels = c("0","1"))
  cm_svm2 = confusionMatrix(svm_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  svm_pred = predict(classifier_svm2, newdata = rad_ext)
  svm_pred = factor(svm_pred, levels = c("0","1"), labels = c("0","1"))
  cm_svm2_ext = confusionMatrix(svm_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 4) NAIVE BAYES CLASSIFIER
  set.seed(42)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  
  classifier_nb = train(form = His ~ .,
                        data = rad_train3,
                        method = "naive_bayes",
                        trControl=trctrl,
                        tuneLength = 20
  )
  
  classifier_nb$bestTune
  classifier_nb
  
  # Predicting the validation set results
  nb_pred = predict(classifier_nb, newdata = rad_val3)
  nb_pred = factor(nb_pred, levels = c("0","1"), labels = c("0","1"))
  cm_nb = confusionMatrix(nb_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  nb_pred = predict(classifier_nb, newdata = rad_ext)
  nb_pred = factor(nb_pred, levels = c("0","1"), labels = c("0","1"))
  cm_nb_ext = confusionMatrix(nb_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 5) Random Forest
  set.seed(42)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  mtry <- sqrt(ncol(rad_train3))
  tunegrid <- expand.grid(mtry = c(2,3,mtry))
  #tunegrid <- expand.grid(mtry = 2)
  classifier_rf = train(form = His ~ .,
                        data = rad_train3,
                        method = "rf",
                        trControl=trctrl,
                        tuneGrid = tunegrid
                        #tuneLength = 20
  )
  
  rf_pred = predict(classifier_rf, newdata = rad_val3[-ncol(rad_val3)])
  rf_pred = factor(rf_pred, levels = c("0","1"), labels = c("0","1"))
  cm_rf = confusionMatrix(rf_pred, reference = rad_val3$His, positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.### 6) ARTIFICIAL NEURAL NETWORK
  
  rf_pred = predict(classifier_rf, newdata = rad_ext[-ncol(rad_ext)])
  rf_pred = factor(rf_pred, levels = c("0","1"), labels = c("0","1"))
  cm_rf_ext = confusionMatrix(rf_pred, reference = rad_ext$His, positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.### 6) ARTIFICIAL NEURAL NETWORK
  
  ### 6) ARTIFICIAL NEURAL NETWORK
  set.seed(42)
  trctrl <- trainControl(method = "repeatedcv", number = 1, repeats = 1, verboseIter = FALSE)
  classifier_nnet = train(form = His ~ .,
                          data = rad_train3,
                          method = "nnet",
                          trace=FALSE,
                          trainControl = trctrl,
                          tuneLength = 20,
                          verboseIter=FALSE
  )
  # Predicting the NN validation set results
  nn_pred = predict(classifier_nnet, newdata = rad_val3[-ncol(rad_val3)])
  nn_pred = factor(nn_pred, levels = c("0","1"), labels = c("0","1"))
  cm_nnet = confusionMatrix(nn_pred, as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  nn_pred = predict(classifier_nnet, newdata = rad_ext[-ncol(rad_ext)])
  cm_nnet_ext = confusionMatrix(nn_pred , reference = rad_ext$His) #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 7) XGBOOST - with grid search
  set.seed(10)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
  trgrid <- expand.grid(nrounds = 1000,
                        eta = c(0.3, 0.1, 0.01, 0.001),
                        max_depth = c(2, 4, 6, 8),
                        gamma = c(0, 1),
                        subsample = c(0.5, 0.6),
                        colsample_bytree=c(0.5,0.6),
                        min_child_weight = seq(1)
  )
  classifier_xgb = train(form = His ~ .,
                         data = rad_train3,
                         method = "xgbTree",
                         trControl = trctrl,
                         #tuneLength = 10
                         tuneGrid = trgrid
  )
  xgb_pred = predict(classifier_xgb, newdata = rad_val3[-ncol(rad_val3)])
  cm_xgb = confusionMatrix(xgb_pred, reference = rad_val3$His) #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  xgb_pred = predict(classifier_xgb, newdata = rad_ext[-ncol(rad_ext)])
  cm_xgb_ext = confusionMatrix(xgb_pred , reference = rad_ext$His) #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 8) PARTIAL LEAST SQUARES - with grid search
  set.seed(42)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  
  classifier_pls = train(form = His ~ .,
                         data = rad_train3,
                         method = "pls",
                         trainControl = trctrl,
                         tuneLength = 20
  )
  # Predicting the validation set results
  pls_pred = predict(classifier_pls, newdata = rad_val3[-ncol(rad_val3)])
  cm_pls = confusionMatrix(pls_pred, reference = rad_val3$His) #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  pls_pred = predict(classifier_pls, newdata = rad_ext[-ncol(rad_ext)])
  cm_pls_ext = confusionMatrix(pls_pred, reference = rad_ext$His) #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 9) Ridge Regression - with grid search
  
  #First make a set of lambda values
  lambda <- 10^seq(-3, 3, length = 100)
  #lambda <- c(0, 0.0001)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  
  set.seed(123)
  classifier_ridge <- train( His ~.,
                             data = rad_train3,
                             method = "glmnet",
                             family = "binomial",
                             trControl = trctrl,
                             tuneGrid = expand.grid(alpha = 0, lambda = lambda)
  )
  
  # Model coefficients
  coef(classifier_ridge$finalModel, classifier_ridge$bestTune$lambda)
  
  classifier_ridge$bestTune
  
  #Make predictions on the validation set
  ridge_pred = predict(classifier_ridge, newdata = rad_val3[-ncol(rad_val3)])
  ridge_pred = factor(ridge_pred, levels = c("0","1"), labels = c("0","1"))
  cm_ridge = confusionMatrix(ridge_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ridge_pred = predict(classifier_ridge, newdata = rad_ext[-ncol(rad_ext)])
  ridge_pred = factor(ridge_pred, levels = c("0","1"), labels = c("0","1"))
  cm_ridge_ext = confusionMatrix(ridge_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation set.
  
  ### 10) LASSO Regression - with grid search
  
  #First make a set of lambda values
  lambda <- 10^seq(-3, 3, length = 100)
  
  set.seed(123)
  classifier_lasso <- train(His ~.,
                            data = rad_train3, 
                            method = "glmnet",
                            family = "binomial",
                            trControl = trainControl("repeatedcv", number = 10, repeats = 3),
                            tuneGrid = expand.grid(alpha = 1, lambda = lambda)
  )
  
  # Make predictions on the validation set
  lasso_pred = predict(classifier_lasso, newdata = rad_val3[-ncol(rad_val3)])
  lasso_pred = factor(lasso_pred, levels = c("0","1"), labels = c("0","1"))
  cm_lasso = confusionMatrix(lasso_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a conf### 11) Elastic Net Regression - with grid search
  
  lasso_pred = predict(classifier_lasso, newdata = rad_ext[-ncol(rad_ext)])
  lasso_pred = factor(lasso_pred, levels = c("0","1"), labels = c("0","1"))
  cm_lasso_ext = confusionMatrix(lasso_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a conf### 11) Elastic Net Regression - with grid search
  
  ##### 11) Elastic Net 
  #First make a set of lambda values
  lambda <- 10^seq(-3, 3, length = 100)
  
  set.seed(123)
  classifier_elastic <- train(His ~.,
                              data = rad_train3, 
                              method = "glmnet",
                              trControl = trainControl("repeatedcv", number = 10, repeats = 3),
                              tunelength = expand.grid(alpha = 0.5, lambda = lambda)
  )
  # Make predictions on the validation set
  elastic_pred = predict(classifier_elastic, newdata = rad_val3[-ncol(rad_val3)])
  elastic_pred = factor(elastic_pred, levels = c("0","1"), labels = c("0","1"))
  cm_elastic = confusionMatrix(elastic_pred, reference = as.factor(rad_val3$His), positive = "1")
  
  elastic_pred = predict(classifier_elastic, newdata = rad_ext[-ncol(rad_ext)])
  elastic_pred = factor(elastic_pred, levels = c("0","1"), labels = c("0","1"))
  cm_elastic_ext = confusionMatrix(elastic_pred, reference = as.factor(rad_ext$His), positive = "1")
  
  #cm_glm, cm_knn2, cm_svm2, cm_nb, cm_rf, cm_nnet, cm_xgb, cm_pls, cm_ridge, cm_lasso, cm_elastic
  
  #tmp <- list(cm_elastic$byClass[7], cm_glm$byClass[7], cm_lasso$byClass[7], cm_ridge$byClass[7], cm_nb$byClass[7], cm_knn2$byClass[7], cm_nnet$byClass[7], cm_rf$byClass[7], cm_svm2$byClass[7], cm_xgb$byClass[7], cm_pls$byClass[7])
  #mp <- list(cm_elastic$overall[1], cm_glm$byClass[1], cm_lasso$byClass[1], cm_ridge$byClass[1], cm_nb$byClass[1], cm_knn2$byClass[7],cm_rf$byClass[7], cm_svm2$byClass[7], cm_pls$byClass[7])
  #tmp <- list(cm_elastic$overall[1], cm_glm$overall[1], cm_lasso$overall[1], cm_ridge$overall[1], cm_knn2$overall[1],cm_rf$overall[1], cm_svm2$overall[1], cm_pls$overall[1])
  
  #tmp <- c(cm_elastic$overall[1], cm_glm$overall[1], cm_lasso$overall[1], cm_ridge$overall[1], cm_knn2$overall[1],cm_rf$overall[1], cm_svm2$overall[1], cm_pls$overall[1])
  
  tmp <- c(cm_elastic$overall[1], cm_glm$overall[1], cm_lasso$overall[1], cm_ridge$overall[1], cm_knn2$overall[1],cm_rf$overall[1], cm_svm2$overall[1], cm_pls$overall[1], cm_xgb$overall[1], cm_nnet$overall[1], cm_nb$overall[1])
  
  #tmp_ext <- list(cm_elastic_ext$overall[1], cm_glm_ext$overall[1], cm_lasso_ext$overall[1], cm_ridge_ext$overall[1], cm_knn2_ext$overall[1],cm_rf_ext$overall[1], cm_svm2_ext$overall[1], cm_pls_ext$overall[1])
  #tmp_ext <- c(cm_elastic_ext$overall[1], cm_glm_ext$overall[1], cm_lasso_ext$overall[1], cm_ridge_ext$overall[1], cm_knn2_ext$overall[1],cm_rf_ext$overall[1], cm_svm2_ext$overall[1], cm_pls_ext$overall[1])

  tmp_ext <- c(cm_elastic_ext$overall[1], cm_glm_ext$overall[1], cm_lasso_ext$overall[1], cm_ridge_ext$overall[1], cm_knn2_ext$overall[1],cm_rf_ext$overall[1], cm_svm2_ext$overall[1], cm_pls_ext$overall[1], cm_xgb_ext$overall[1], cm_nnet_ext$overall[1], cm_nb_ext$overall[1])
  
  res <- cbind(res,tmp)
  res_ext <- cbind(res_ext,tmp_ext)
  
}

## model validation results as a DF 
res <- as.data.frame(res)

### Rank features in order of importance according to the one of the models
importance <- varImp(classifier_nb, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance, top=33, main="Feature Importance for NB model - predicting Mallignancy", ylab="Features")

breaksList = seq(0.5, 1, by = 0.01)

rownames(res) <- c("ENet", "GLM", "LASSO","Ridge", "KNN", "RF", "SVM", "PLS", "XGB", "NNet", "NB")
colnames(res) <- c("Pearson", "Spearman", "Kendall", "MI", "LASSO", "ENet", "Boruta", "RFE","PCA")

## heatmap of validation results- will need to add colnames/ rownames to make work
pheatmap(res,display_numbers = TRUE, treeheight_row = 0, treeheight_col =0, fontsize = 15, order=FALSE,  cluster_rows=F, cluster_cols=F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)

rownames(res_ext) <- c("ENet", "GLM", "LASSO","Ridge", "KNN", "RF", "SVM", "PLS", "XGB", "NNet","NB")
colnames(res_ext) <- c("Pearson", "Spearman", "Kendall", "MI", "LASSO", "ENet", "Boruta", "RFE","PCA")

pheatmap(res_ext,display_numbers = TRUE, treeheight_row = 0, treeheight_col =0, order=FALSE,  cluster_rows=F, cluster_cols=F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)

##################################Best Stats Method###################################
# 
# #Select Boruta x KNN
# i = 7
# 
# rad_train3 <- NULL
# rad_train3 <- harm_rad_train2[ ,unlist(features_select[i])]
# rad_val3 <- NULL
# rad_val3 <- harm_rad_val2[ ,unlist(features_select[i])]
# rad_ext <- NULL
# rad_ext <- extdata[ ,unlist(features_select[i])]
# 
# rad_train3$His <- y
# rad_val3$His <- val_set_full$His 
# synth_train <- rad_train3
# rad_ext$His <- extdata$His
# 
# #synth_train <- ADASYN(rad_train3, perc_min = 50, k=5)
# #synth_val <- ADASYN(rad_val3, perc_min = 50, k=5)
# 
# #rad_train3 <- synth_train
# ###rad_val3 <- synth_val
# 
# ###Fitting logistic regression to the training set
# rad_train3$His <- as.factor(rad_train3$His)
# rad_val3$His <- as.factor(rad_val3$His)
# rad_ext$His <- as.factor(rad_ext$His)
# 
# ## 2 KNN
# trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# set.seed(42)
# 
# ###Fitting knn to the training set
# classifier_knn2 = train(form = His ~ .,
#                         data = rad_train3,
#                         method = "knn",
#                         trControl=trctrl,
#                         tuneLength = 20)
# 
# classifier_knn2$bestTune #Predicting the validation set results
# knn2_pred = predict(classifier_knn2, newdata = rad_val3[-ncol(rad_val3)])
# knn2_pred = factor(knn2_pred, levels = c("0","1"), labels = c("0","1"))
# cm_knn2 = confusionMatrix(knn2_pred, reference = as.factor(rad_val3$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation
# 
# knn2_pred = predict(classifier_knn2, newdata = rad_ext[-ncol(rad_ext)])
# knn2_pred = factor(knn2_pred, levels = c("0","1"), labels = c("0","1"))
# cm_knn2_ext = confusionMatrix(knn2_pred, reference = as.factor(rad_ext$His), positive = "1") #Create a confusion matrix comparing the predicted results to the real outcome in validation


##LASSO

i = 2

rad_train4 <- NULL
rad_train4 <- harm_rad_train2[ ,unlist(features_select[i])]
rad_val4 <- NULL
rad_val4 <- harm_rad_val2[ ,unlist(features_select[i])]
rad_ext4 <- NULL
rad_ext4 <- extdata[ ,unlist(features_select[i])]

#rad_train4$His <- y
#rad_val4$His <- val_set_full$His 
#synth_train4 <- rad_train4
#rad_ext4$His <- extdata$EGFRMut

nFolds <- 50
foldid <- sample(rep(seq(nFolds), length.out = nrow(rad_train4)))

forlasso <- rad_train4

list.of.fits <- list()
for (i in 0:10){
  cvfit.name <- paste0("cvalpha", i/10)
  #list.of.fits[[cvfit.name]]<- cv.glmnet(x=as.matrix(forlasso), y=rad_train4$His, type.measure = "mse", alpha = i/10, family = "binomial", nfolds = nFolds, foldid = foldid)
  list.of.fits[[cvfit.name]]<- cv.glmnet(x=as.matrix(forlasso), y=train_set_full$His, type.measure = "mse", alpha = i/10, family = "binomial", nfolds = nFolds, foldid = foldid)
}

results <- data.frame()

#test it
for (i in 0:10){
  fit.name <- paste0("alpha", i/10)
  cvfit.name <- paste0("cvalpha", i/10)
  list.of.fits[[fit.name]] <- glmnet(x=as.matrix(forlasso), y=train_set_full$His, type.measure = "mse", alpha=i/10,  family="binomial", nfolds = nFolds, lambda=list.of.fits[[cvfit.name]]$lambda.min, foldid = foldid)
  
  predicted <- predict(list.of.fits[[fit.name]], newx=as.matrix(rad_val4, list.of.fits[[cvfit.name]]$lambda.min))
  
  mse <- mean((y= val_set_full$His - predicted)^2)
  results <- rbind(results, mse)
}

cvfit <- cv.glmnet(x=as.matrix(forlasso), y=train_set_full$His, type.measure = "mse", alpha = 1, family = "binomial", nfolds = nFolds, foldid = foldid)
plot(cvfit)

#Get cross validated R squared for goodness of fit
rsq = 1 - cvfit$cvm/var(y)
plot(cvfit$lambda,rsq)

fit <- glmnet(x=as.matrix(forlasso), y=train_set_full$His, type.measure = "mse", alpha = 0.1, family = "binomial", nfolds = nFolds, lambda=cvfit$lambda.min, foldid = foldid)
fit$beta[,1]

coef(fit)
tidy(fit)

#Prediction
prediction_model_dis <- predict(fit, newx=as.matrix(rad_train4, cvfit$lambda.min))

#internal validation
prediction_model_val <- predict(fit, newx=as.matrix(rad_val4, cvfit$lambda.min))

#external validation and testing sets
prediction_ext <- predict(fit, newx=as.matrix(rad_ext4, cvfit$lambda.min))

#EV_nor<-cbind(EV_nor,prediction_EV)
#T_nor<-cbind(T_nor,prediction_T)
#M_nor<-cbind(M_nor,prediction_M)

##################################AUROC Plots###################################

library(ROCit)

rad_train4$His <- train_set_full$His 
rad_val4$His <- val_set_full$His 
synth_train4 <- rad_train4
rad_ext4$His <- extdata$His

cate_dis <- as.numeric(rad_ext4$His)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

score_dis <- prediction_ext

ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("EGFR-RPV: EGFR Exclusive Positivity - External Testing")

cate_dis <- as.numeric(rad_train4$His)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

score_dis <- prediction_model_dis

ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("EGFR-RPV: EGFR Exclusive Positivity - Training")

##################################Survival Plots################################

############prediction k means function######################
predict.kmeans <- function(object, newdata){
  centers <- object$centers
  n_centers <- nrow(centers)
  dist_mat <- as.matrix(dist(rbind(centers, newdata)))
  dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
  max.col(-dist_mat)
}

########################clustering - Ordinary Lasso#################
suppressPackageStartupMessages(install_github("michaelway/ggkm"))
#gitcreds::gitcreds_delete()
install.packages("survminer")
library(ggkm)
library(survival)
library(survminer)

#Survival data available
T_final <- cbind(extdata, prediction_ext)
#T_final <- T_final[complete.cases(T_final[, 8]),]

EV_final <- cbind(val_set_full, prediction_model_val)
#EV_final <- EV_final[complete.cases(EV_final[, 3]),]

Dis_final <- cbind(train_set_full, prediction_model_dis)
#Dis_final <- Dis_final[complete.cases(Dis_final[, 3]),]

Dis_combined <- merge(Dis_final, EV_final, all=TRUE)
rad_combined <- merge(rad_train4, rad_val4, all=TRUE)

no_partition <- 5

#Needs to rerun lines 883-934
####Combined Discovery Cohort######################
lassomat <- rad_combined
sig_array <- fit$beta@i
postlasso <- lassomat[, fit$beta@i]
postlasso <- cbind(postlasso, Dis_combined$s0)

set.seed(20)
km.res <- kmeans(postlasso, no_partition, nstart = 1)
Dis_combined$RPV_Grouping <- km.res$cluster

####Training Cohort######################
lassomat <- rad_train4
sig_array <- fit$beta@i
postlasso <- lassomat[, fit$beta@i]
postlasso <- cbind(postlasso, prediction_model_dis)

set.seed(12)
km.res <- kmeans(postlasso, no_partition, nstart = 1)
Dis_final$RPV_Grouping <- km.res$cluster

#Interal Validation
lassomat <- rad_val4
postlasso <- lassomat[, fit$beta@i]
postlasso <- cbind(postlasso, prediction_model_val)

set.seed(11)
km.res <- kmeans(postlasso, no_partition, nstart = 1)
EV_final$RPV_Grouping <- km.res$cluster

#External Testing
lassomat <- rad_ext4
postlasso <- lassomat[, fit$beta@i]
postlasso <- cbind(postlasso, prediction_ext)

set.seed(100)
no_partition <- 2
km.res <- kmeans(postlasso, no_partition, nstart = 1)
T_final$RPV_Grouping <- km.res$cluster


#EV_final$LCRPV_Grouping <- predict(km.res, postlasso_val)

Dis_final <- Dis_combined

ggkm(survfit(Surv(Overall.survival..days.,OS.event)~RPV_Grouping,data= Dis_final),main = "Kaplan-Meier Plots of Groups Stratified Using \n LC-RPV - Discovery Cohort", margins=c(12,15), pval=T,table=T, legend = T, legendposition = c(0.85, 0.15), ystrataname = "Group", xlab = "Time-to-Event (Days)", cex = 0.8, marks=T)

###Only needed if no. clusters > 2
#re-sort groups

#6 partition, 3,4 as one
for (i in 1:nrow(Dis_final)){
  #if ((Dis_final$RPV_Grouping[i] == 1) | (Dis_final$RPV_Grouping[i] == 3)) {
  if (Dis_final$RPV_Grouping[i] == 2) {
    Dis_final$RPV_Grouping[i] <- 1
  }
  else{
    Dis_final$RPV_Grouping[i] <- 2
  }
}

Dis_final$`EGFR-RPV`<-Dis_final$RPV_Grouping

coxph(Surv(Overall.survival..days.,OS.event) ~ `EGFR-RPV`, data= Dis_final) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

fit <- survfit(Surv(Overall.survival..days.,OS.event) ~ RPV_Grouping, data= Dis_final)

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()) 

#Fancy KM Plot
a<-ggsurvplot(
  fit,     # survfit object with calculated statistics.
  data = Dis_final,               # data used to fit survival curves. 
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1095),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations
  # in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  title = "Kaplan-Meier Plots of Stratified Groups By EGFR-RPV",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000, 1095))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 200, 400, 600, 800, 1000),]) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 200, 400, 600, 800, 1000),]) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

##############Nomogram Creation#######################
#install.packages(c("PASWR","rms"))
library(PASWR)
library("rms")

T_roc <- Dis_combined[ -c(10,11,12) ]
T_roc <- subset(T_roc, T_roc$Histology<3)

T_roc <- T_roc[complete.cases(T_roc$His),]

T_roc$`EGFR-RPV` <- abs(T_roc$s0)
t.data <- datadist(T_roc)
options(datadist = 't.data')

fit <- lrm(formula = His ~ Gender + Ethnicity + `EGFR-RPV`, data = T_roc)
plot(nomogram(fit, fun = function(x)plogis(x)))

#################Draw Forest Plot for Multivariable Regression##################
#install.packages("xfun")
#install.packages("forestplot")
#install.packages("forestmodel")
library(forestplot)
library(caret)

T_roc$Gender <- ifelse(T_roc$Gender==1,"Male","Female")
T_roc$Ethnicity<-ifelse(T_roc$Ethnicity == 1,"Caucasian",ifelse(T_roc$Ethnicity == 3,"Asian","Other"))
#T_roc$Histology<-ifelse(T_roc$Histology == 1,"SqC", "Adeno")

set.seed(100)

logitmod1 <- glm(formula = His ~ Gender, family = binomial(link = "logit"), data = T_roc)

summary(logitmod)
confint(logitmod)

#+as.factor(Ethnicity)+as.factor(Histology)+LCRPV

library(tidyverse)
library(forestmodel)

plot1 <- forest_model(logitmod1)
plot1

logitmod2 <- glm(formula = His ~ Ethnicity, family = binomial(link = "logit"), data = T_roc)
#logitmod3 <- glm(formula = His ~ Histology, family = binomial(link = "logit"), data = T_roc)
logitmod4 <- glm(formula = His ~ `EGFR-RPV`, family = binomial(link = "logit"), data = T_roc)
plot2 <- forest_model(logitmod2)
#plot3 <- forest_model(logitmod3)
plot4 <- forest_model(logitmod4)

logitmod <- glm(formula = His ~ Gender + Ethnicity + `EGFR-RPV`, family = binomial(link = "logit"), data = T_roc)


library(ggpubr)

ggarrange(plot1, plot2, plot4, ncol=1, nrow=3)
