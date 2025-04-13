library("mpath")
library("zic")
library("caret")
library("glmnet")
library("pscl")
library("forecast")
library("splines")
library("stringr")

data("docvisits")
# PUT FOLLOWING GRAPH IN REPORT:
# barplot(with(docvisits, table(docvisits)), ylab = "Frequency", xlab = "Doctor office visits")

rename_polynomials <- function(degree, name) {
    # no intercept so degree = num of variables
    names <- numeric(degree)
    for (i in 1:degree) {
        names[i] <- paste(name, i, sep = "_")
    }
    return(names)
}

n <- nrow(docvisits)
age <- bs(docvisits$age, 3)[1:n, ]
colnames(age) <- rename_polynomials(3, "age")
hlth <- bs(docvisits$health, 3)[1:n, ]
colnames(hlth) <- rename_polynomials(3, "hlth")
hdeg <- bs(docvisits$hdegree, 3)[1:n, ]
colnames(hdeg) <- rename_polynomials(3, "hdeg")
schl <- bs(docvisits$schooling, 3)[1:n, ]
colnames(schl) <- rename_polynomials(3, "schl")
hhin <- bs(docvisits$hhincome, 3)[1:n, ]
colnames(hhin) <- rename_polynomials(3, "hhin")

zipdoc <- cbind.data.frame(docvisits$docvisits, age, hlth, hdeg, schl, hhin, handicap, married, children, self, civil, bluec, employed, public, addon)
names(zipdoc)[1] <- "docvisits"
yvar <- names(zipdoc)[1]
xvars <- names(zipdoc)[-1]
zvars <- xvars

# same groups as Chatterjee et al.
est_groups <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

# make this a function that takes in the gooogle or gooogleplus function as a parameter.
real_data_run <- function(method, n) {
    mae <- 0
    mase <- 0

    for (partition_seed in 1:n) {
        print(partition_seed)
        set.seed(partition_seed)
        train_indices <- createDataPartition(zipdoc$docvisits, p = 0.8, list = FALSE)
        zipdoc_train <- zipdoc[train_indices, ]
        zipdoc_test <- zipdoc[-train_indices, ]
        print(system.time(fit <- method(data = zipdoc_train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = "grLasso", dist = "poisson", group = est_groups)))

        betahat <- fit$coefficients$count
        gammahat <- fit$coefficients$zero

        # Calculate predicted phi
        z.test <- as.matrix(cbind(1, zipdoc_test[, zvars]))
        phi.hat <- 1 / (1 + exp(-z.test %*% gammahat))  # Calculate phi.hat

        # Calculate predicted lambda
        x.test <- as.matrix(cbind(1, zipdoc_test[, xvars]))
        lam.hat <- exp(x.test %*% betahat)

        # Calculate predicted y
        y.pred <- (1 - phi.hat) * lam.hat

        # Extract actual y values
        y.test <- zipdoc_test[, yvar]
        y.train <- zipdoc_train[, yvar]

        # Create forecast object
        forecast <- structure(list(mean = y.pred, fitted = y.test, x = y.train), class = "forecast")

        # Calculate accuracy measures (MCC and AUC) and round to 4 decimals
        measures <- c(round(accuracy(forecast, y.test)[2, c(3, 6)], 4))
        mae <- mae + measures[1]
        mase <- mase + measures[2]
    }

    return(list(
        mae = mae / n,
        mase = mase / n
    ))
}

gooogle_result <- real_data_run(gooogle, 10)
gooogleplus_result <- real_data_run(gooogleplus, 10)
print(gooogle_result)
print(gooogleplus_result)

test_gp_time <- function() {
    train_indices <- createDataPartition(zipdoc$docvisits, p = 0.8, list = FALSE)
    zipdoc_train <- zipdoc[train_indices, ]
    print(system.time(gooogleplus(data = zipdoc_train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = "grLasso", dist = "poisson", group = est_groups)))
}

test_gp_time()

res <- zeroinfl(docvisits ~ ., data = zipdoc, dist = "negbin")
reszip <- zeroinfl(docvisits ~ ., data = zipdoc, dist = "poisson")
print(reszip$coefficients)
names(res)
names(reszip$theta)
