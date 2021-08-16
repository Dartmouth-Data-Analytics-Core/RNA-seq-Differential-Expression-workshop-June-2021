## Day2-03 Linear Modeling

# read in the example data
dat <- read.csv("data/lm-example-data.csv", stringsAsFactors=FALSE)

# explore it quickly
head(dat)
str(dat)

# plot
plot(dat$gene_exp ~ dat$hba1c,
     ylab = "Expression (Gene X)",
     xlab = "Hba1c score",
     main = "Gene X exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# fit a linear model with gene expression as the response
lm1 <- lm(dat$gene_exp ~ dat$hba1c)
lm1


# generate plot again
plot(dat$gene_exp ~ dat$hba1c,
     ylab = "Expression (Gene X)",
     xlab = "Hba1c score",
     main = "Gene X exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# add the model on the scatterplot
abline(lm1, lty=2)

# calculate the predicted gene expression values using the model
pre <- predict(lm1)

# plot the difference between the predicted and the true values
segments(dat$hba1c, dat$gene_exp, dat$hba1c, pre,
         col="cornflowerblue")
#### Note: These are the residuals!


### hypothesis testing with linear models

sum_lm1 <- summary(lm1)
sum_lm1

# get the coefficients table
coef(sum_lm1)

# get the coefficients themselves
coef(sum_lm1)[,1]

# get the P-value for the hba1c coefficient
coef(sum_lm1)[2,4]


## Gene Y example

# read in the example data
dat2 <- read.csv("data/lm-example-data-geneY.csv", stringsAsFactors=FALSE)

# plot
plot(dat2$gene_exp ~ dat2$hba1c,
     ylab = "Expression (Gene Y)",
     xlab = "Hba1c score",
     main = "Gene Y exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# fit a linear model with gene expression as the response
lm1 <- lm(dat2$gene_exp ~ dat2$hba1c)
summary(lm1)
pre <- predict(lm1)

# add the model on the scatterplot
abline(lm1, lty=2)

# plot the difference between the predicted and the true values
segments(dat2$hba1c, dat2$gene_exp, dat2$hba1c, pre, col="cornflowerblue")


## Categorical variables and linear modeling

# read in the example data
dat3 <- read.csv("data/lm-example-3.csv", stringsAsFactors=FALSE, row.names = 1)

# quickly explore it
head(dat3)
table(dat3$subject_group)
# Note: Controls are coded as 0, cases are coded as 1

# visualize the data
boxplot(dat3$exp_geneX ~ dat3$subject_group ,
        ylab = "Expression (Gene X)",
        xlab = "Subject group",
        main = "Gene X exp. in subject groups",
        col = c("indianred", "cornflowerblue"), pch = 16, las = 1)


# run the linear model and evaluate
lm_2 <- lm(dat3$exp_geneX ~ dat3$subject_group)
summary(lm_2)
