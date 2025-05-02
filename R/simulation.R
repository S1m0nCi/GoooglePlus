library(rockchalk)
library(mpath)
library(outliers)
library(forecast)
library(gamlss.dist)
library(utils)

# Constants for Simulation
START_SEED <- 1000
NUM_SIMULATIONS <- 200
NUM_TRAIN <- c(200, 500, 1000)
NUM_TEST <- c(80, 200, 400)
RHO <- 0.1 # correlation-related parameter in [0, 1]
PHI <- 0.4 # controls zero inflation component intercept term in (0, 1)
BETA <- c(-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2, 8), rep(0, 24))
GAMMA <- c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2, 8), rep(0, 24))
GROUP_STRUCT <- rep.int(8, 5) # element i is the number of variables in the ith group
GOOOGLE_PENALTIES <- c("grLasso", "grMCP", "grSCAD")
GOOOGLEPLUS_PENALTIES <- c("grLasso", "grMCP", "grSCAD", "grALasso")

SAMEGRP.OVERLAP <- TRUE

# Derived Constants
GROUPS <- numeric(sum(GROUP_STRUCT))
counter <- 0
for (i in seq_along(GROUP_STRUCT)) {
    for (j in 1:(GROUP_STRUCT[i])) {
        GROUPS[counter + j] <- i
    }
    counter <- counter + GROUP_STRUCT[i]
}

# Simulation Functions
gen_zip_data <- function(
    n.train,
    n.test,
    grpsize,
    rho,
    phi,
    seedval,
    beta,
    gamma
) {
    # Define beta
    beta <- c(5, beta)

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
    counter <- 0
    for (g in 1:ngrp) {
        for (j in 1:grpsize[g]) {
            X[, counter + j] <- (Z[, g] + R[, counter + j]) / sqrt(2)
        }
        counter <- counter + grpsize[g]
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
    sensspec_count <- sens_spec(betahat, beta)
    sensspec_zero <- sens_spec(gammahat, gamma)
    sens_count <- sensspec_count$sensitivity
    spec_count <- sensspec_count$specificity
    sens_zero <- sensspec_zero$sensitivity
    spec_zero <- sensspec_zero$specificity

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

    # Calculate accuracy measures (MCC and AUC) and round to 4 decimal places
    measures <- c(round(accuracy(forecast, y.test)[2, c(3, 6)], 4),
                  sensitivity = round(sens, 4),
                  specificity = round(spec, 4),
                  sens_count = round(sens_count, 3),
                  spec_count = round(spec_count, 3),
                  sens_zero = round(sens_zero, 3),
                  spec_zero = round(spec_zero, 3))

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
        time.taken <- system.time(fit.summary <- fit.method(data = train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = method, dist = "poisson", group = group, samegrp.overlap = SAMEGRP.OVERLAP))

        # Predict measures for the test set
        predict.measures <- measures.func(train = train, test = test, fit = fit.summary, yvar = yvar, xvars = xvars, zvars = zvars, beta = count, gamma = zero)

        # Append measures and time to the matrix
        measures.mat <- rbind(measures.mat, c(predict.measures, time.taken[3]))
    }

    # Calculate standard errors using b.mean with bootstrapping
    measures.se <- t(apply(apply(measures.mat[, c(1, 2)], 2, function(x) return(as.numeric(x))), 2, b.mean, num = 1000, na.rm = TRUE))

    # Calculate medians for mae and mase, not sens and spec
    accuracy.median <- apply(measures.mat[, c(1:2, 9)], 2, function(x) { median(x, na.rm = TRUE) })
    vs.mean <- apply(measures.mat[, 3:8], 2, function(x) { mean(x) })

    # Summarize MAE with median and standard error
    mae.summary <- paste(round(accuracy.median[1], 3), "(", round(measures.se[1], 3), ")", sep = "")

    # Summarize MASE with median and standard error
    mase.summary <- paste(round(accuracy.median[2], 3), "(", round(measures.se[2], 3), ")", sep = "")

    # Combine results into output structure
    output <- c(MAE = mae.summary, MASE = mase.summary, vs.mean[1], vs.mean[2], time.taken = round(accuracy.median[3], 2), vs.mean[3], vs.mean[4], vs.mean[5], vs.mean[6])
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

get_groups <- function(group_struct) {
    groups <- numeric(sum(group_struct))
    counter <- 0
    for (i in seq_along(group_struct)) {
        for (j in 1:(group_struct[i])) {
            groups[counter + j] <- i
        }
        counter <- counter + group_struct[i]
    }
    return(groups)
}

run_sim <- function(sim_set) {
    data_list <- lapply(1:sim_set$num_sim, function(i) {
        return(gen_zip_data(sim_set$num_train, sim_set$num_test, sim_set$group_struct, sim_set$rho, sim_set$phi, START_SEED + i, sim_set$beta, sim_set$gamma))
    })
    groups <- get_groups(sim_set$group_struct)
    results_df <- data.frame(
        id = integer(),
        package = character(),
        MAE = numeric(),
        MASE = numeric(),
        sensitivity = numeric(),
        specificity = numeric(),
        time = numeric(),
        rho = numeric(),
        phi = numeric(),
        num_sim = integer(),
        num_train = integer(),
        num_test = integer(),
        sens_count = numeric(),
        spec_count = numeric(),
        sens_zero = numeric(),
        spec_zero = numeric()
    )
    for (penalty in GOOOGLEPLUS_PENALTIES) {
        if (penalty %in% GOOOGLE_PENALTIES) {
            g_measures <- measures.summary(gooogle, sim_set$num_train, data_list, penalty, groups)
            g_df <- data.frame(
                id = sim_set$id,
                package = "gooogle",
                penalty = penalty,
                MAE = g_measures[1],
                MASE = g_measures[2],
                sensitivity = g_measures[3],
                specificity = g_measures[4],
                time = g_measures[5],
                rho = sim_set$rho,
                phi = sim_set$phi,
                num_sim = sim_set$num_sim,
                num_train = sim_set$num_train,
                num_test = sim_set$num_test,
                sens_count = g_measures[6],
                spec_count = g_measures[7],
                sens_zero = g_measures[8],
                spec_zero = g_measures[9]
            )
            results_df <- rbind(results_df, g_df)
        }
        gp_measures <- measures.summary(gooogleplus, sim_set$num_train, data_list, penalty, groups)
        gp_df <- data.frame(
            id = sim_set$id,
            package = "gooogleplus",
            penalty = penalty,
            MAE = gp_measures[1],
            MASE = gp_measures[2],
            sensitivity = gp_measures[3],
            specificity = gp_measures[4],
            time = gp_measures[5],
            rho = sim_set$rho,
            phi = sim_set$phi,
            num_sim = sim_set$num_sim,
            num_train = sim_set$num_train,
            num_test = sim_set$num_test,
            sens_count = gp_measures[6],
            spec_count = gp_measures[7],
            sens_zero = gp_measures[8],
            spec_zero = gp_measures[9]
        )
        results_df <- rbind(results_df, gp_df)
    }
    return(results_df)
}

run_simulations <- function(set_list, write_file) {
    sim_results <- data.frame(
        id = integer(),
        package = character(),
        MAE = numeric(),
        MASE = numeric(),
        sensitivity = numeric(),
        specificity = numeric(),
        time = numeric(),
        rho = numeric(),
        phi = numeric(),
        num_sim = integer(),
        num_train = integer(),
        num_test = integer(),
        sens_count = numeric(),
        spec_count = numeric(),
        sens_zero = numeric(),
        spec_zero = numeric()
    )
    for (i in seq_along(set_list)) {
        print(set_list[[i]]$id)
        sim_set <- set_list[[i]]
        sim_result <- run_sim(sim_set)
        sim_results <- rbind(sim_results, sim_result)
    }
    print(paste("C:\\Users\\gbsim\\Documents\\R-scripts\\GoooglePlus\\R\\", write_file, ".csv", sep = ""))
    write.table(sim_results, file = paste("C:\\Users\\gbsim\\Documents\\R-scripts\\GoooglePlus\\R\\", write_file, ".csv", sep = ""), row.names = FALSE, sep = ",")
}

run_large_simulation <- function(large_sim_set, write_file) {
    sim_results <- data.frame(
        id = integer(),
        package = character(),
        MAE = numeric(),
        MASE = numeric(),
        sensitivity = numeric(),
        specificity = numeric(),
        time = numeric(),
        rho = numeric(),
        phi = numeric(),
        num_sim = integer(),
        num_train = integer(),
        num_test = integer(),
        sens_count = numeric(),
        spec_count = numeric(),
        sens_zero = numeric(),
        spec_zero = numeric()
    )
    for (i in seq_along(large_sim_set)) {
        print(large_sim_set[[i]]$id)
        sim_set <- large_sim_set[[i]]
        for (i in seq_along(NUM_TRAIN)) {
            sim_set$num_train <- NUM_TRAIN[i]
            sim_set$num_test <- NUM_TEST[i]
            sim_result <- run_sim(sim_set)
            sim_results <- rbind(sim_results, sim_result)
        }
    }
    write.table(sim_results, file = paste("C:\\Users\\gbsim\\Documents\\R-scripts\\GoooglePlus\\R\\", write_file, ".csv", sep = ""), row.names = FALSE, sep = ",")
}


initial_set <- list(
    list(
        id = 1, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0, 0, 0, 0),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0, 0, 0, 0),
        group_struct = c(1, 2, 1, 2, 2)
    ),
    list(
        id = 2, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0, 0, 0, 0),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, 0, 0, 0, 0),
        group_struct = c(1, 2, 1, 3, 2, 2)
    ),
    list(
        id = 3, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4)
    ),
    list(
        id = 4, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4)
    ),
    list(
        id = 5, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4)
    ),
    list(
        id = 6, # to identify
        num_sim = 2,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    )
)

refining_simulation_set <- list(
    list(
        id = 11, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 12, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.5,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 13, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.9,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 21, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 22, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.5,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 23, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.9,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 31, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.1,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 32, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.5,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 33, # to identify
        num_sim = 20,
        num_train = 200,
        num_test = 50,
        rho = 0.9,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    )
)

large_simulation_set <- list(
    list(
        id = 11, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.1,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 12, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.5,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 13, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.9,
        phi = 0.3,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 21, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.1,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 22, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.5,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 23, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.9,
        phi = 0.4,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 31, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.1,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 32, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.5,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    ),
    list(
        id = 33, # to identify
        num_sim = NUM_SIMULATIONS,
        num_train = NUM_TRAIN,
        num_test = NUM_TEST,
        rho = 0.9,
        phi = 0.5,
        # below here can be identified with id: not to include in report tables, talk about outside
        beta = c(-1, -0.5, -0.25, -0.1, 0.3, 0.4, -0.4, 0.1, -0.3, 0.2, -0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        gamma = c(-0.4, -0.3, -0.2, -0.1, 0.2, -0.6, 0.2, -0.2, -0.3, 0.4, 0.5, rep(0, 4), rep(0, 4), rep(0, 4)),
        group_struct = c(1, 2, 1, 3, 4, 4, 4, 4)
    )
)

# run_simulations(initial_set, "tp")

# run_large_simulation(large_simulation_set, "results_sgt")

# SAMEGRP.OVERLAP <- FALSE

# run_large_simulation(large_simulation_set, "results_sgf")

# list(
#     id = 1, # to identify
#     num_sim = 10,
#     num_train = 200,
#     num_test = 50,
#     rho = 0.1,
#     phi = 0.4,
#     # below here can be identified with id: not to include in report tables, talk about outside
#     beta = c(-1, -0.5, -0.25, -0.1),
#     gamma = c(-0.4, -0.3, -0.2, -0.1),
#     group_struct = c(4)
# ),

# list(
#     id = 6, # to identify
#     num_sim = 10,
#     num_train = 200,
#     num_test = 50,
#     rho = 0.1,
#     phi = 0.4,
#     # below here can be identified with id: not to include in report tables, talk about outside
#     beta = c(-1, -0.5, -0.25, -0.1),
#     gamma = c(-0.4, -0.3, -0.2, -0.1),
#     group_struct = c(1, 2, 1)
# )
