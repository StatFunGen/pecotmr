set.seed(985115)
n <- 350
p <- 16
sigmasq_error <- 0.5
zeroes <- rbinom(p, 1, 0.6)
beta.true <- rnorm(p, 1, sd = 4)
beta.true[zeroes] <- 0

X <- cbind(matrix(rnorm(n * p), nrow = n))
X <- scale(X, center = TRUE, scale = FALSE)
y <- X %*% matrix(beta.true, ncol = 1) + rnorm(n, 0, sqrt(sigmasq_error))
y <- scale(y, center = TRUE, scale = FALSE)

# Calculate sufficient statistics
XtX <- t(X) %*% X
Xty <- t(X) %*% y
yty <- t(y) %*% y

# Set the prior
K <- 9
sigma0 <- c(0.001, .1, .5, 1, 5, 10, 20, 30, .005)
omega0 <- rep(1 / K, K)

# Calculate summary statistics
b.hat <- sapply(1:p, function(j) {
    summary(lm(y ~ X[, j]))$coefficients[-1, 1]
})
s.hat <- sapply(1:p, function(j) {
    summary(lm(y ~ X[, j]))$coefficients[-1, 2]
})
R.hat <- cor(X)
var_y <- var(y)
sigmasq_init <- 1.5

# Run PRS CS
maf = rep(0.5, length(b.hat)) # fake MAF
LD <- list(blk1 = R.hat)
sumstats <- list(BETA = b.hat, MAF = maf)
write.table(
    data.frame(sumstats),
    "sumstats.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)
write.table(LD$blk1, "LD.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

out2 <- pecotmr::prs_cs(b.hat, LD, n, maf = maf, verbose = TRUE, seed = 999)
print(out2$beta_est)

#Rscript generate_test_for_prs_cs.R
#Running Markov Chain Monte Carlo (MCMC) sampler...
#Iteration  100 of 1000
#Iteration  200 of 1000
#Iteration  300 of 1000
#Iteration  400 of 1000
#Iteration  500 of 1000
#Iteration  600 of 1000
#Iteration  700 of 1000
#Iteration  800 of 1000
#Iteration  900 of 1000
#Iteration 1000 of 1000
#Estimated global shrinkage parameter: 0.573608
#MCMC sampling completed.
#            [,1]
# [1,]  0.2014556
# [2,]  1.0877955
# [3,] -1.5719068
# [4,]  1.1523724
# [5,] -0.8181124
# [6,] -0.8719565
# [7,] -1.5151863
# [8,]  1.7940113
# [9,]  2.0403245
#[10,] -1.7775449
#[11,] -0.8175015
#[12,]  6.1956311
#[13,] -0.6421237
#[14,] -1.9703685
#[15,]  0.8190806
#[16,] -0.3304925
