library("mpath")
library("zic")
library("caret")
library("glmnet")
library("pscl")
library("forecast")


data(docvisits)
# PUT FOLLOWING GRAPH IN REPORT:
# barplot(with(docvisits, table(docvisits)), ylab = "Frequency", xlab = "Doctor office visits")

dt <- docvisits[, -(2:3)]
tmp <- model.matrix(~ age30 * health + age35 * health + age40 * health + age45 * health + age50 * health + age55 * health + age60 * health, data = dt)[, -(1:9)]
dat <- cbind(dt, tmp)
zipdoc <- dat[, -(15:21)]

yvar <- names(zipdoc)[1]
xvars <- names(zipdoc)[-1]
zvars <- xvars

est_groups <- c(1, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4, 1, 1, 5, 5, 5, 5, 5, 5, 5)
# group 1: health, handicap, addon, public (insurance details)
# group 2: hdegree, schooling, hhincome (social status)
# group 3: married, children (personal life)
# group 4: self, civil, bluec, employed (work type)
# group 5: age30TRUE:health, health:age35TRUE, health:age40TRUE, health:age45TRUE, health:age50TRUE, health:age55TRUE, health:age60TRUE (age)
# gooogle

zipdoc_train_df <- as.data.frame(zipdoc_train)

# make this a function that takes in the gooogle or gooogleplus function as a parameter.
real_data_run <- function(method, n) {
    mae <- 0
    mase <- 0

    for (partition_seed in 1:n) {
        set.seed(partition_seed)
        train_indices <- createDataPartition(zipdoc$docvisits, p = 0.8, list = FALSE)
        zipdoc_train <- zipdoc[train_indices, ]
        zipdoc_test <- zipdoc[-train_indices, ]
        fit <- method(data = zipdoc_train_df, yvar = yvar, xvars = xvars, zvars = zvars, penalty = "grLasso", dist = "poisson", group = est_groups)

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

gooogle_result <- real_data_run(gooogle, 100)
gooogleplus_result <- real_data_run(goooglePlus, 100)
print(gooogle_result)
print(gooogleplus_result)
