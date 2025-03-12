library(rockchalk)
library(mpath)
library(outliers)
library(forecast)
library(gamlss.dist)

gen_zip_data <- function(
    n.train,
    n.test,
    grpsize,
    rho,
    phi,
    seedval
) {
    # Define beta
    beta <- c(-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2, 8), rep(0, 24))
    beta <- c(5, beta)

    # Define gamma
    gamma <- c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2, 8), rep(0, 24))

    # Adjust gamma based on phi
    if (phi == 0.3) gamma <- c(-1, gamma)
    if (phi == 0.4) gamma <- c(-0.5, gamma)
    if (phi == 0.5) gamma <- c(0, gamma)

    # Set seed if provided
    if (!is.null(seedval)) set.seed(seedval)

    # Define n, p, and ngrp
    n <- n.train + n.test
    p <- sum(grpsize)
    ngrp <- length(grpsize)

    # Generate random normal matrix (R)
    R <- matrix(rnorm(n * p), n, p)

    # Create correlation matrix (V)
    V <- matrix(0, ngrp, ngrp)
    for (i in 1:ngrp) {
        for (j in 1:ngrp) {
            V[i, j] <- rho^(abs(i - j))
        }
    }

    # Generate multivariate normal matrix (Z)
    Z <- mvrnorm(n, mu = rep(0, ngrp), Sigma = V)

    # Create design matrix (X)
    X <- matrix(0, n, p)
    for (g in 1:ngrp) {
        for (j in 1:grpsize[g]) {
            X[, (g - 1) * grpsize[g] + j] <- (Z[, g] + R[, (g - 1) * grpsize[g] + j]) / sqrt(2)
        }
    }

    # Standardize X
    X <- scale(X)

    # Set column names for X
    colnames(X) <- paste("X", c(1:ncol(X)), sep = "")

    # Define variable names
    xvars <- colnames(X)
    zvars <- xvars

    # Capture zero inflation (if seedval is null)
    if (is.null(seedval)) {
        rn <- round(runif(1) * 10^5)
        set.seed(rn)
        y <- rzi(n = n, x = X, z = X, a = beta, b = gamma, family = "poisson")

        set.seed(rn)
        # Capture from console output
        yout <- capture.output(rzi(n, x = X, z = X, a = beta, b = gamma, family = "poisson"))
        zeroinfl <- as.numeric(substring(yout[1], 15))
    } else {
        y <- rzi(n = n, x = X, z = X, a = beta, b = gamma, family = "poisson")

        # Capture from console output
        yout <- capture.output(rzi(n, x = X, z = X, a = beta, b = gamma, family = "poisson"))
        zeroinfl <- as.numeric(substring(yout[1], 15))
    }

    # Combine data into a data frame
    data <- cbind.data.frame(y, X)

    # Return list
    return(list(data = data,
                yvar = "y",
                xvars = xvars,
                zvars = zvars,
                zeroinfl = zeroinfl,
                coefficients = list(
                    count = beta,
                    zero = gamma
                )))
}

# let us run one simulation
# 40 variables: this is decided in the data generation
output <- gen_zip_data(200, 50, rep.int(8, 5), 0.1, 0.4, 200)
data <- output$data
yvar <- output$yvar
xvars <- output$xvars
zvars <- output$zvars

print(system.time(sim_result_grLasso <- goooglePlus(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), penalty = "grLasso", lambda_min = 0)))
print(system.time(sim_result_grLasso_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grLasso", lambda.min = 0)))
sim_result_grLasso
sim_result_grLasso$coefficients
sim_result_grLasso_old$coefficients

print(sim_result_grLasso_old$coefficients)
print(sim_result_grLasso$params[[100]])
#sim_result_grLasso$params
#sim_result_grMCP <- goooglePlus(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), penalty = "grMCP")
#sim_result_grSCAD <- goooglePlus(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), penalty = "grSCAD")
#print("LASSO")
#print(sim_result_grLasso$coefficients)
#print(sim_result_grLasso$bic)
#print("MCP")
#print(sim_result_grMCP$coefficients)
#print(sim_result_grMCP$bic)
#print("SCAD")
#print(sim_result_grSCAD$coefficients)
#print(sim_result_grSCAD$bic)
#print(sim_result_grLasso_old$coefficients)
#sim_result_grMCP_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grMCP")
#sim_result_grSCAD_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grSCAD")

b.mean <- function(vec, num, na.rm = TRUE) {
    # Remove missing values if na.rm is TRUE
    if (na.rm) vec <- vec[!is.na(vec)]

    # Remove outliers from the vector
    vec <- rm.outlier(vec)
    # Generate resamples using lapply
    resamples <- lapply(1:num, function(i) sample(vec, replace = TRUE))

    # Calculate mean of each resample using sapply
    r.mean <- sapply(resamples, mean)

    # Calculate standard error
    std.err <- sqrt(var(r.mean))

    # Round standard error and return
    return(round(std.err, 2))
}

#' Title
#'
#' @description For a given training and test dataset and the fitted coefficients this function calculates MAE and MASE
#' @param train Training dataset obtained from the original dataset
#' @param test Complement of the training dataset used for prediction
#' @param fit The output of the function fit.method
#' @param yvar Name of the outcome variable
#' @param xvars Name of the predictor variables for the count model
#' @param zvars Name of the predictor variables for the zero model
#'
#' @return The predictive measures MAE and MASE calculated from the function accuracy of the forecast package
measures.func <- function(train, test, fit, yvar, xvars, zvars, beta, gamma) {
    # Check if fit is missing (NA)
    # Extract coefficients
    betahat <- fit$coefficients$count
    gammahat <- fit$coefficients$zero

    # Calculate sensitivity and specificity
    sensspec <- sens_spec(c(betahat, gammahat), c(beta, gamma))
    sens <- sensspec$sensitivity
    spec <- sensspec$specificity

    # Calculate predicted phi
    z.test <- as.matrix(cbind(1, test[, zvars]))
    phi.hat <- 1 / (1 + exp(-z.test %*% gammahat))  # Calculate phi.hat

    # Calculate predicted lambda
    x.test <- as.matrix(cbind(1, test[, xvars]))
    lam.hat <- exp(x.test %*% betahat)

    # Calculate predicted y
    y.pred <- (1 - phi.hat) * lam.hat

    # Extract actual y values
    y.test <- test[, yvar]
    y.train <- train[, yvar]

    # Create forecast object
    forecast <- structure(list(mean = y.pred, fitted = y.test, x = y.train), class = "forecast")

    # Calculate accuracy measures (MCC and AUC) and round to 4 decimals
    measures <- c(round(accuracy(forecast, y.test)[2, c(3, 6)], 4),
                  sensitivity = round(sens, 4),
                  specificity = round(spec, 4))

    # Return the measures
    return(measures)
}

#' Title
#'
#' @description This function fits the ZINB model with different penalties to the training part of a dataset and calculates MAE and MASE from the test set. It outputs the median of MAE and MASE over all the simulated datasets
#' @param fit.method the function being used to fit the model
#' @param n.train Sample size in the training dataset
#' @param data.list output of datagen.sim.all
#' @param method Different penalties
#' @param group Vector containing grouping structure of the covariates
#'
#' @return The median of MAE and MASE, calculated over all the simulated datasets
measures.summary <- function(fit.method, n.train, data.list, method, group) {
    # Suppress warnings during iteration
    options(warn = -1)

    # Initialize measures matrix
    measures.mat <- NULL

    # Iterate over repetitions
    for (i in 1:length(data.list)) {
        # Extract data for current iteration
        dataset <- data.list[[i]]
        train <- dataset$data[1:n.train, ]
        test <- dataset$data[-(1:n.train), ]

        # Extract variables
        yvar <- dataset$yvar
        xvars <- dataset$xvars
        zvars <- dataset$zvars

        # Extract actual coefficients
        count <- dataset$coefficients$count
        zero <- dataset$coefficients$zero

        # Fit the model and capture time
        time.taken <- system.time(fit.summary <- fit.method(data = train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = method, dist = "poisson", group = group))

        # Predict measures for the test set
        predict.measures <- measures.func(train = train, test = test, fit = fit.summary, yvar = yvar, xvars = xvars, zvars = zvars, beta = count, gamma = zero)

        # Append measures and time to the matrix
        measures.mat <- rbind(measures.mat, c(predict.measures, time.taken[3]))
    }
    print(measures.mat)

    # Calculate standard errors using b.mean with bootstrapping
    measures.se <- t(apply(apply(measures.mat[, c(1, 2)], 2, function(x) return(as.numeric(x))), 2, b.mean, num = 1000, na.rm = TRUE))

    # Calculate medians for each measure
    measures.median <- apply(measures.mat, 2, function(x) { median(x, na.rm = TRUE) })

    # Summarize MAE with median and standard error
    mae.summary <- paste(round(measures.median[1], 3), "(", round(measures.se[1], 3), ")", sep = "")

    # Summarize MASE with median and standard error
    mase.summary <- paste(round(measures.median[2], 3), "(", round(measures.se[2], 3), ")", sep = "")

    # Combine results into output structure
    output <- c(MAE = mae.summary, MASE = mase.summary, measures.median[3], measures.median[4], time.taken = measures.median[5])

    # Restore warning settings
    options(warn = 0)

    # Return the summary output
    return(output)
}

sens_spec <- function(estimate, actual) {
    correct_nonzero <- 0
    correct_zero <- 0
    for (i in 1:length(estimate)) {
        if (estimate[i] == 0  && actual[i] == 0) {
            correct_zero <- correct_zero + 1
        }
        else if (estimate[i] != 0 && actual[i] != 0) {
            correct_nonzero <- correct_nonzero + 1
        }
    }
    return(list(
        sensitivity = correct_nonzero / sum(actual != 0),
        specificity = correct_zero / sum(actual == 0)
    )
    )
}

data.list <- lapply(1:100, function(i) {
    return(gen_zip_data(200, 50, rep.int(8, 5), 0.1, 0.4, i))
})

print("Group LASSO")
simulation_results_gooogleplus <- measures.summary(goooglePlus, 200, data.list, "grLasso", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(simulation_results_gooogleplus)
simulation_results_gooogle <- measures.summary(gooogle, 200, data.list, "grLasso", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(simulation_results_gooogle)

print("Group MCP")
gooogleplus_mcp <- measures.summary(goooglePlus, 200, data.list, "grMCP", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(gooogleplus_mcp)
gooogle_mcp <- measures.summary(gooogle, 200, data.list, "grMCP", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(gooogle_mcp)

print("Group SCAD")
gooogleplus_scad <- measures.summary(goooglePlus, 200, data.list, "grSCAD", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(gooogleplus_scad)
gooogle_scad <- measures.summary(gooogle, 200, data.list, "grSCAD", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(gooogle_scad)

print("Group ALASSO")
gooogleplus_alasso <- measures.summary(goooglePlus, 200, data.list, "grALasso", c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(gooogleplus_alasso)
