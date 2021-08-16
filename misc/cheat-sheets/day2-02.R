# Day2-02 Statistics for DE

# generate an empty vector to store P-values in
p.value <- c()

# generate 2 random variables (r.v.s) 1000 times and run a t.test on each pair of r.v.s
# the r.v.s will have the same mean and distribution
for(i in 1:1000){
  # simulate random variables
  x <- rnorm(n = 20, mean = 0, sd = 1)
  y <- rnorm(n = 20, mean = 0, sd = 1)
  # run t.test and extract p-value
  p.value[i] <- t.test(x, y)$p.value
}

# count number of P-values less than our alpha
table(p.value < 0.05)

# order vector of P-values
p.value <- p.value[order(p.value)]

# visualize P-value magnitude
plot(-log10(p.value), las = 1,
     col = "cornflowerblue",
     main = "-log 10 P-values",
     xlab = "ordered P-values")
#### Note: taking -log10 of small number gives big number!

# add a horizontal line
abline(h=-log10(0.05), col = "red", lty = 2)

# add some useful labels
text(600, 1.5, "Small P-values higher up")
text(600, 1.1, "Large P-values lower down")

### Bonferroni correction
# visualize P-value magnitude
plot(-log10(p.value), las = 1,
     col = "cornflowerblue",
     main = "-log 10 P-values",
     xlab = "ordered P-values", ylim = c(0,5))

# add a horizontal line
abline(h=-log10(0.05), col = "red", lty = 2)
text(600, 4.5, "Bonferroni")

# add a horizontal line
abline(h=-log10(0.05/1000), col = "red", lty = 2)
text(600, 1.5, "Original threshold")

## adjust P-vals with Bonferroni corrections

# bonferroni correction
p.adj.bonf <- p.adjust(p.value, method = "bonferroni")
p.adj.bonf

# check if any are signifciant
table(p.adj.bonf < 0.05)

### FDR corrections

BiocManager::install("qvalue")
library(qvalue)
qvalue(p.value)$qvalues
p.adj.fdr <- qvalue(p.value)$qvalues

# check how many are sig.
table(p.adj.fdr < 0.05)
