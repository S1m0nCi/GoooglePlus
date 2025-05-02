library("mpath")
library("zic")
library("caret")
library("glmnet")
library("pscl")
library("forecast")
library("splines")
library("stringr")

GOOOGLE_PENALTIES <- c("grLasso", "grMCP", "grSCAD")
GOOOGLEPLUS_PENALTIES <- c("grLasso", "grMCP", "grSCAD", "grALasso")
FILE_PATH <- "C:\\Users\\gbsim\\Documents\\R-scripts\\GoooglePlus\\R\\"

data(docvisits)
# PUT FOLLOWING GRAPH IN REPORT:
# barplot(with(docvisits, table(docvisits)), ylab = "Frequency", xlab = "Doctor office visits", cex.lab = 1.5, cex.axis = 1.5)

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

zipdoc <- cbind.data.frame(docvisits$docvisits, age, hlth, hdeg, schl, hhin,
                           docvisits$handicap, docvisits$married, docvisits$children,
                           docvisits$self, docvisits$civil, docvisits$bluec,
                           docvisits$employed, docvisits$public, docvisits$addon)
names(zipdoc)[1] <- "docvisits"
yvar <- names(zipdoc)[1]
xvars <- names(zipdoc)[-1]
zvars <- xvars

# same groups as Chatterjee et al.
groups <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

# make this a function that takes in the gooogle or gooogleplus function as a parameter.
real_data_run <- function(method, penalty, n) {
    mae <- numeric(5 * n)
    mase <- numeric(5 * n)
    time <- 0
    coefficients_list <- vector("list", length = n * 5)
    print(as.character(substitute(method)))
    print(penalty)

    for (partition_seed in 1:n) {
        print(partition_seed)
        set.seed(partition_seed)
        folds <- createFolds(zipdoc$docvisits, k = 5)
        for (i in 1:5) {
            zipdoc_train <- zipdoc[-folds[[i]], ]
            zipdoc_test <- zipdoc[folds[[i]], ]
            run_time <- system.time(fit <- method(data = zipdoc_train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = penalty,
                                                  dist = "poisson", group = groups))
            time <- time + run_time[3]
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

            mae[partition_seed * 5 + i] <- measures[1]
            mase[partition_seed * 5 + i] <- measures[2]
            coefficients_list[[partition_seed * 5 + i]] <- fit$coefficients
        }
    }

    mae_order <- order(mae)
    median_coefficients <- coefficients_list[[median(mae_order)]]

    return(list(
        mae = median(mae),
        mase = median(mase),
        time = time / (n * 5),
        coefficients = median_coefficients
    ))
}

run_all <- function(method, penalties, n, write_file) {
    results <- data.frame(penalty = character(),
                          mae = numeric(),
                          mase = numeric(),
                          time = numeric())
    coefficients_df <- data.frame(row.names = 1:(length(groups) + 1))

    for (penalty in penalties) {
        run_results <- real_data_run(method, penalty, n)
        pen_results <- data.frame(penalty = penalty,
                                  mae = run_results$mae,
                                  mase = run_results$mase,
                                  time = run_results$time)
        results <- rbind.data.frame(results, pen_results)
        coefficients_df <- cbind.data.frame(coefficients_df, run_results$coefficients$count)
        coefficients_df <- cbind.data.frame(coefficients_df, run_results$coefficients$zero)
    }
    write.table(results, file = paste(FILE_PATH, write_file, "_results", ".csv", sep = ""), row.names = FALSE, sep = ",")
    write.table(coefficients_df, file = paste(FILE_PATH, write_file, "_coeff", ".csv", sep = ""), row.names = FALSE, sep = ",")
}

# run_all(gooogle, GOOOGLE_PENALTIES, 100, "real_data_gooogle")
# run_all(gooogleplus, GOOOGLEPLUS_PENALTIES, 100, "real_data_gooogleplus")
# gooogle_result <- real_data_run(gooogle, "grLasso", 10)
# gooogleplus_result <- real_data_run(gooogleplus, "grLasso", 10)
# gooogle_result
# gooogleplus_result
