library(pscl)
library(vscDebugger)
library(glmnet)

goooglePlus <-  function(
    data,
    xvars,
    zvars,
    yvar,
    group = 1:ncol(data),
    samegrp.overlap = TRUE,
    penalty = c("grLasso", "grMCP", "grSCAD"),
    dist = c("poisson"),
    nlambda = 100,
    lambda,
    lambda_min = ifelse((nrow(data[, unique(c(xvars, zvars))]) > ncol(data[, unique(c(xvars, zvars))])), 1e-4, .05),
    lambda_max,
    log.lambda = TRUE,
    crit = "BIC",
    alpha = 1,
    eps = .001,
    max.iter = 1000,
    gmax = length(unique(group)),
    gamma = ifelse(penalty == "grSCAD", 4, 3),
    warn = TRUE
)
{
    # getting the appropriate columns from the data frame - 1d output only
    y <- data.frame(data[, yvar])
    X <- data.frame(data[, xvars])
    # getting all of the explanatory variables: count-influencing and zero-inflation-influencing
    pred.names <- union(xvars, zvars)
    # renaming
    names(X) <- paste(names(X), ".count", sep = "")
    names(y) <- yvar

    # remove zero from groups so that we can set the intercept group equal to 0
    if (0 %in% group) {
        group_min <- min(group)
        group <- group + rep(group_min + 1, length(group))
    }

    # gets the groups of each explanatory variable into a list of group nums format
    group.x <- group[which(pred.names %in% xvars)]  # group of xvars

    if (is.null(zvars)) {  # if there is no covariate in the zero model
        Z <- NULL
        group.z <- NULL
        data <- cbind.data.frame(y, X)
        xvars <- names(X)
        fit.formula <- as.formula(paste(yvar, "~", paste(paste(xvars, collapse = "+"), "|1"), sep = ""))
    } else {
        Z <- data.frame(data[, zvars])
        # rename Z dataframe
        names(Z) <- paste(names(Z), ".zero", sep = "")

        # the below code just ensures no overlap of X and Z group numbers: they are treated as completely separate groups: hence we add the max of groups x to ensure none are shared
        if (samegrp.overlap) {  # if X and Z assign same groups for shared covariates
            group.z <- group[which(pred.names %in% zvars)]
        } else {
            group.z <- max(group.x) + group[which(pred.names %in% zvars)]
        }
        # cbind - bind colums (rbind - bind rows)
        data <- cbind.data.frame(y, X, Z)
        xvars <- names(X)
        zvars <- names(Z)
        # "|" separates count model | zero inflation model
        fit.formula <- as.formula(paste(yvar, "~", paste(paste(xvars, collapse = "+"), "|", paste(zvars, collapse = "+")), sep = ""))
    }

    p <- length(xvars)
    q <- length(zvars)

    # vcov matrix also includes for the intercepts
    fit.zero <- zeroinfl(fit.formula, dist = dist, data = data)

    # keep the intercepts in the mle and the final data, which we want to have p + q + 2 columns [Zeng]
    # different to Gooogle where the approximation is used.
    b2.mle <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)

    fit.coefficients <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)
    zeroinf_residuals <- fit.zero$residuals
    # get the covariance matrix needed for transforming the data
    vcov <- fit.zero$vcov

    names(fit.coefficients)[1] <- "count.intercept"
    names(fit.coefficients)[p + 2] <- "zero.intercept"
    # get the pseudo data factor vcov (sigma) to the -1/2
    evcov <- eigen(vcov)
    # vcov is symmetric square matrix therefore the matrix of its eigenvectors should be orthogonal ie AT = A^-1
    inv_sqrt_diag <- diag(vapply(evcov$values, function(num) {
        num ^ (-1 / 2)
    }, numeric(1)), nrow = nrow(evcov$vectors), ncol = ncol(evcov$vectors))
    pseudo_factor <- evcov$vectors %*% inv_sqrt_diag %*% t(evcov$vectors)

    pseudo_X <- pseudo_factor
    pseudo_Y <- as.vector(pseudo_factor %*% b2.mle)

    # reorder the group index and the covariates
    group <- c(0, group.x, 0, group.z)

    # initialise pseudo_X reordered by group
    group_ordered_indices <- order(group)
    # get the indices for putting into unordered form again (original form)
    unordering_indices <- numeric(length(group))
    for (i in 1:length(group_ordered_indices)) {
        unordering_indices[group_ordered_indices[i]] <- i
    }
    pseudo_X <- pseudo_X[, group_ordered_indices, drop = FALSE]
    fit.coefficients <- fit.coefficients[group_ordered_indices]

    group <- sort(group)
    unique_groups <- unique(group)
    # scaling and group sorting has now already been done
    # now we implement the group coordinate descent algorithm
    # according to Zeng, we do not need to transform back to non-LSA, so non-pseudo data
    # however we will need to orthogonalise and then transform back

    # get the group matrices
    # use 'match'
    group_start_indices <- c(match(unique_groups, group), length(group) + 1)
    penalty_factors <- numeric(length(unique_groups))
    # the intercepts are to be unpenalised - setting penalisation to 0.
    penalty_factors[1] <- 0

    # we want penalty factors of zero for the intercepts, so the last two elements of the group
    for (i in 2:(length(group_start_indices) - 1)) {
        # group_matrices[[i]] <- pseudo_X[, group_start_indices[i]:group_start_indices[i + 1] - 1]
        penalty_factors[i] <- sqrt(group_start_indices[i + 1] - group_start_indices[i])
    }
    # group_matrices[length(unique_groups)] <- pseudo_X[, group_start_indices[length(unique_groups)]:ncol(pseudo_X)]

    # orthonormalise the group matrices
    group_matrices_orth_full <- orthonormalise_full_matrix(pseudo_X, group_start_indices)
    group_matrix_orthonormal <- group_matrices_orth_full[[1]]
    orthonormalisation_factors <- group_matrices_orth_full[[2]]

    # Multiple lambda values - one parameter only, which is lambda
    if (missing(lambda)) {
        if (missing(lambda_max)) {
            # we may need glmnet here instead of zeroinfl because our method relies on group coordinate descent. Just for finding lambda_max
            group_matrix_df <- as.data.frame(group_matrix_orthonormal)
            names(group_matrix_df) <- names(fit.coefficients)
            intercepts_fit <- glm(
                as.formula(paste(yvar, "~ 0 + count.intercept + zero.intercept")),
                data = cbind(data.frame(y = unlist(pseudo_Y)), group_matrix_df)
            )
            # will give matrix of 2 cols, instead of vector
            # make a matrix of X, Z and rearrange for groups as necessary
            zmax <- maxgrad(group_matrix_orthonormal, intercepts_fit$residuals, group_start_indices, penalty_factors) / nrow(group_matrix_orthonormal)
            lambda_max <- zmax
        }
        if (log.lambda) {
            # creating in descending order
            if (lambda_min == 0) {
                lambda <- c(exp(seq(log(lambda_max), log(0.001 * lambda_max), length = nlambda - 1)), 0)
            } else {
                lambda <- exp(seq(log(lambda_max), log(lambda_min * lambda_max), length = nlambda))
            }
        }
        else {
            if (lambda_min == 0) {
                lambda <- c(seq(lambda_max, 0.001 * lambda_max, length = nlambda - 1), 0)
            } else {
                lambda <- seq(lambda_max, lambda_min * lambda_max, length = nlambda)
            }
        }
    }

    orthonormalised_coefficients_list <- group_coordinate_descent(
        group_matrix_orthonormal,
        lambda,
        penalty_factors,
        fit.coefficients,
        zeroinf_residuals,
        eps,
        group_start_indices,
        pseudo_Y,
        penalty = penalty
    )

    # Step 1: De-orthonormalise and unorder coefficients
    coefficients_list <- lapply(orthonormalised_coefficients_list, function(coeff) {
        return(deorthonormalise_coefficients(
            coeff,
            orthonormalisation_factors,
            group_start_indices
        )[unordering_indices])
    })

    mat_X <- as.matrix(cbind(1, X))
    mat_Z <- as.matrix(cbind(1, Z))
    mat_y <- as.matrix(y)

    # Step 2: Find bic and optimise
    coefficient_log_likelihood <- vapply(coefficients_list, function(coeff) {
        return(
            zip_log_likelihood(
                mat_X,
                mat_Z,
                coeff[1:(p + 1)],
                coeff[(p + 2): (p + q + 2)],
                mat_y
            )
        )
    }, numeric(1))

    # get number of model parameters (non-zero coefficients) and exclude intercepts - we have two intercepts
    # at the moment, this may be overstated => Zeng's approach of e^lambda may be better.
    df <- vapply(coefficients_list, function(x) { return(sum(x != 0) - 2) }, numeric(1))
    coefficients_bic <- vapply(1:length(coefficient_log_likelihood), function(i) {
        return(df[i] * log(nrow(pseudo_X)) - 2 * coefficient_log_likelihood[i])
    }, numeric(1))

    min_bic_index <- which.min(coefficients_bic)
    optimal_coefficients <- coefficients_list[[min_bic_index]]
    optimal_lambda <- lambda[i]
    optimal_bic <- coefficients_bic[min_bic_index]
    count_coefficients <- optimal_coefficients[1:(p + 1)]
    zero_coefficients <- optimal_coefficients[(p + 2):(p + q + 2)]
    names(count_coefficients) <- names(fit.zero$coefficients$count)
    names(zero_coefficients) <- names(fit.zero$coefficients$zero)
    opt_params <- list(count = count_coefficients, zero = zero_coefficients)

    return(list(
                params = coefficients_list,
                group = group,
                lambda = lambda,
                df = df,
                loss = coefficient_log_likelihood,
                bic = coefficients_bic,
                penalty = penalty,
                n = nrow(X),
                # iter = iterations,
                coefficients = opt_params,
                opt_lambda = optimal_lambda,
                opt_bic = optimal_bic))
}

orthonormalise_full_matrix <- function(all_groups_matrix, group_start_indices) {
    num_cols <- ncol(all_groups_matrix)
    orthonormalised_matrix <- matrix(, nrow = nrow(all_groups_matrix), ncol = num_cols)
    matrix_factors <- vector("list", length(group_start_indices))
    for (i in 1:(length(group_start_indices) - 1)) {
        # get the matrix:
        group_matrix <- all_groups_matrix[, group_start_indices[i]:(group_start_indices[i + 1] - 1), drop = FALSE]
        orth_group_matrix_and_factor <- orthonormalise_group_matrix(group_matrix)
        orthonormalised_group_matrix <- orth_group_matrix_and_factor[[1]]
        matrix_factors[[i]] <- orth_group_matrix_and_factor[[2]]
        orthonormalised_matrix[, group_start_indices[i]:(group_start_indices[i + 1] - 1)] <- orthonormalised_group_matrix
    }
    return(list(orthonormalised_matrix, matrix_factors))
}

orthonormalise_group_matrix <- function(matrix) {
    gram <- (1 / nrow(matrix)) * t(matrix) %*% matrix
    # gram is a square matrix
    eigenresult <- eigen(gram, symmetric = TRUE)
    eigenvector_matrix <- eigenresult$vectors
    inv_diag_matrix <- diag(vapply(eigenresult$values, function(num) {
        num ^ (-1 / 2)
    }, numeric(1)), nrow = nrow(gram), ncol = ncol(gram))
    matrix_factor <- eigenvector_matrix %*% inv_diag_matrix
    # keep the matrix factor for transformation back later on.
    matrix_and_factor <- vector("list", 2)
    matrix_and_factor[[1]] <- matrix %*% matrix_factor
    matrix_and_factor[[2]] <- matrix_factor
    return(matrix_and_factor)
}

deorthonormalise_coefficients <- function(coeff, matrix_factors, group_start_indices) {
    deorth_coeff <- numeric(length(coeff))
    for (i in 1:(length(group_start_indices) - 1)) {
        deorth_coeff[group_start_indices[i]:(group_start_indices[i + 1] - 1)] <- matrix_factors[[i]] %*% coeff[group_start_indices[i]:(group_start_indices[i + 1] - 1)]
    }
    return(deorth_coeff)
}

l2_distance <- function(u, v) {
    return(sqrt(sum((u - v) ^ 2)))
}

group_coordinate_descent <- function(
    group_mat_orth, # orthonormal matrix of all groups/variables
    lambda,
    penalty_factors,
    initial_b,
    initial_residuals,
    eps,
    group_start_indices,
    y,
    max_iterations = 1000,
    gamma = ifelse(penalty == "grSCAD", 4, 3),
    penalty = c("grLasso", "grMCP", "grSCAD")
) {
    # we need to output all of the coefficients from the different lambdae
    # we then calculate loss etc using bic later
    # also need to worry about how to reorder
    # 'coefficients' below may need to be a vector instead of a list
    coefficients_list <- vector("list", length(lambda))
    n <- nrow(group_mat_orth)
    penalty <- match.arg(penalty)

    if (penalty == "grLasso") {
        threshold <- vector_soft_threshold
    } else if (penalty == "grMCP") {
        threshold <- vector_mcp_firm_threshold
    } else if (penalty == "grSCAD") {
        threshold <- vector_scad_firm_threshold
    }

    for (i in 1:length(lambda)) {
        lam <- lambda[i]

        # do not keep coefficients in group structure
        standard_deviation_y <- sqrt(sum(y ** 2) / nrow(group_mat_orth))
        convergence_threshold <- eps * standard_deviation_y
        penalty_factors <- penalty_factors * lam
        r <- y - group_mat_orth %*% initial_b
        current_b <- initial_b
        for (iterations in 1:max_iterations) {
            # do the intercepts first. We treat them as separate, not as one group. We just 'store' them as one group.
            max_change <- 0
            lambda_0 <- penalty_factors[1]
            for (j in group_start_indices[1]:(group_start_indices[2] - 1)) {
                b <- current_b[j]
                X <- group_mat_orth[, group_start_indices[j], drop = FALSE] # will not be changed
                z <- (1 / n) * t(X) %*% r + b
                new_b <- threshold(z, lambda_0, gamma)
                b_change <- new_b - b
                r <- r - X %*% (b_change) # is there an error in the paper by Breheny, which says to use t(X)?
                b <- new_b
                current_b[j] <- b
                max_change <- max(max_change, abs(b_change))
            }

            for (j in 2:(length(group_start_indices) - 1)) {
                lambda_i <- penalty_factors[j]
                b <- current_b[group_start_indices[j]:(group_start_indices[j + 1] - 1)]
                X <- group_mat_orth[, group_start_indices[j]:(group_start_indices[j + 1] - 1), drop = FALSE] # will not be changed
                z <- (1 / n) * t(X) %*% r + b
                new_b <- threshold(z, lambda_i, gamma)
                b_change <- new_b - b
                r <- r - X %*% (b_change) # is there an error in the paper by Breheny, which says to use t(X)?
                b <- new_b
                current_b[group_start_indices[j]:(group_start_indices[j + 1] - 1)] <- b
                max_change <- max(max_change, max(abs(b_change)))
            }

            if (max_change < convergence_threshold) {
                break
            }
        }
        coefficients_list[[i]] <- current_b
    }
    return(coefficients_list)
}

soft_threshold <- function(num, lambda) {
    if (num > lambda) {
        return(num - lambda)
    } else if (num < - lambda) {
        return(lambda - num)
    }
    return(0)
}

vector_soft_threshold <- function(vect, lambda, gamma) {
    norm <- sqrt(sum(vect^2))
    return(soft_threshold(norm, lambda) * vect / norm)
}

mcp_firm_threshold <- function(num, lambda, gamma) {
    if (abs(num) <= lambda * gamma) {
        return(soft_threshold(num, lambda) / (1 - 1 / gamma))
    } else {
        return(num)
    }
}

scad_firm_threshold <- function(num, lambda, gamma) {
    if (abs(num) <= 2 * lambda) {
        return(soft_threshold(num, lambda))
    } else if (abs(num) <= gamma * lambda) {
        return(soft_threshold(num, lambda * gamma / (gamma - 1)) / (1 - 1 / (gamma - 1)))
    } else {
        return(num)
    }
}

vector_mcp_firm_threshold <- function(vect, lambda, gamma) {
    norm <- sqrt(sum(vect^2))
    return(mcp_firm_threshold(norm, lambda, gamma) * vect / norm)
}

vector_scad_firm_threshold <- function(vect, lambda, gamma) {
    norm <- sqrt(sum(vect^2))
    return(scad_firm_threshold(norm, lambda, gamma) * vect / norm)
}

maxgrad <- function(X, residuals, group_start_indices, penalty_factors) {
    # XX_orth, r, K1, group_multipliers
    zmax <- 0

    # ignore the intercepts here as they are not affected by lambda, as they are not penalised
    for (i in 2:(length(group_start_indices) - 1)) {
        # get group size
        group_length <- group_start_indices[i + 1] - group_start_indices[i]
        Z <- numeric(group_length)

        # + 1 because we are using j - start
        for (j in group_start_indices[i]:(group_start_indices[i + 1] - 1)) {
            # dot product of column and residuals
            # ie how we get XjT * residuals vector, but we are getting each one by one
            Z[j + 1 - group_start_indices[i]] <- sum(X[, j] * residuals)
        }
        # Z is equal to XjT * r
        # ^ assuming that coefficient is zero

        # now, when we penalise we are using group size-derived penalty factors, lambda_j = lambda * sqrt(group_size)
        # this is what actually gets used in the coordinate descent, so we have to take account of that
        # sqrt(sum(Z^2)) is the norm of Z
        # so we need sqrt(sum(Z^2)) < lambda_j (we have to use for group Xj) => sqrt(sum(Z^2)) < lambda * sqrt(group_size)[or penalty factor]
        # hence we get the below lower bound for lambda to penalise (and converge!) to zero. We want to find the greatest lower bound that we can.
        z <- sqrt(sum(Z^2)) / (penalty_factors[i])
        zmax <- max(zmax, z)
    }
    return(zmax)
}

mean_squared_error <- function(pseudo_X, coeff, actual_obs) {
    predictions <- pseudo_X %*% coeff
    return(sum((predictions - actual_obs)^2) / nrow(pseudo_X))
}

# Taken from Lambert's paper
zip_log_likelihood <- function(X, Z, x_coeff, z_coeff, y) {
    # indices of where y is zero
    indices_y0 <- which(y == 0)
    indices_yg0 <- which(y > 0)
    X_0 <- X[indices_y0, ]
    X_g0 <- X[indices_yg0, ]
    Z_0 <- Z[indices_y0, ]
    y_g0 <- y[indices_yg0, ]

    sum_log_y0 <- sum(log(exp(Z_0 %*% z_coeff) + exp(-exp(X_0 %*% x_coeff))))

    sum_all_y <- sum(log(1 + exp(Z %*% z_coeff)))

    sum_yg0 <- sum(y_g0 * (X_g0 %*% x_coeff) - exp(X_g0 %*% x_coeff))

    # sum over all y > 0 of (log(y!) = sum to y from 1 of logn)
    # we could possibly make this even faster
    log_sum <- cumsum(log(1:max(y_g0)))
    sum_yg0_factorial <- sum(log_sum[y_g0])

    return(sum_log_y0 + sum_yg0 - sum_all_y - sum_yg0_factorial)
}

library(mpath)
library(rockchalk)

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
    return(list(data = data, yvar = "y", xvars = xvars, zvars = zvars, zeroinfl = zeroinfl))
}

# let us run one simulation
# 40 variables: this is decided in the data generation
output <- gen_zip_data(200, 50, rep.int(8, 5), 0.1, 0.4, 200)
data <- output$data
yvar <- output$yvar
xvars <- output$xvars
zvars <- output$zvars

#print(system.time(sim_result_grLasso <- goooglePlus(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), penalty = "grLasso")))
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
library(Gooogle)
#print(system.time(sim_result_grLasso_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grLasso")))
#print(sim_result_grLasso_old$coefficients)
#sim_result_grMCP_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grMCP")
#sim_result_grSCAD_old <- gooogle(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), dist = "poisson", penalty = "grSCAD")
library(outliers)
library(forecast)

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
measures.func <- function(train, test, fit, yvar, xvars, zvars) {
    # Check if fit is missing (NA)
    # Extract coefficients
    betahat <- fit$coefficients$count
    gammahat <- fit$coefficients$zero

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
    measures <- c(round(accuracy(forecast, y.test)[2, c(3, 6)], 4))


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
#' @param ITER Number of simulations
#' @param group Vector containing grouping structure of the covariates
#'
#' @return The median of MAE and MASE, calculated over all the simulated datasets
measures.summary <- function(fit.method, n.train, data.list, method, ITER, group) {
    # Suppress warnings during iteration
    options(warn = -1)

    # Initialize measures matrix
    measures.mat <- NULL

    # Iterate over repetitions
    for (i in 1:ITER) {
        # Extract data for current iteration
        dataset <- data.list[[i]]
        train <- dataset$data[1:n.train, ]
        test <- dataset$data[-(1:n.train), ]

        # Extract variables
        yvar <- dataset$yvar
        xvars <- dataset$xvars
        zvars <- dataset$zvars

        # Fit the model and capture time
        time.taken <- system.time(fit.summary <- fit.method(data = train, yvar = yvar, xvars = xvars, zvars = zvars, penalty = method, dist = "poisson", group = group))

        # Predict measures for the test set
        predict.measures <- measures.func(train = train, test = test, fit = fit.summary, yvar = yvar, xvars = xvars, zvars = zvars)

        # Append measures and time to the matrix
        measures.mat <- rbind(measures.mat, c(predict.measures, time.taken[3]))
    }

    # Calculate standard errors using b.mean with bootstrapping
    print(measures.mat)
    measures.se <- t(apply(apply(measures.mat[, -3], 2, function(x) return(as.numeric(x))), 2, b.mean, num = 1000, na.rm = TRUE))

    # Calculate medians for each measure
    measures.median <- apply(measures.mat, 2, function(x) { median(x, na.rm = TRUE) })

    # Summarize MAE with median and standard error
    mae.summary <- paste(round(measures.median[1], 3), "(", round(measures.se[1], 3), ")", sep = "")

    # Summarize MASE with median and standard error
    mase.summary <- paste(round(measures.median[2], 3), "(", round(measures.se[2], 3), ")", sep = "")

    # Combine results into output structure
    output <- c(MAE = mae.summary, MASE = mase.summary, time.taken = measures.median[3])

    # Restore warning settings
    options(warn = 0)

    # Return the summary output
    return(output)
}

data.list <- lapply(1:10, function(i) {
    return(gen_zip_data(200, 50, rep.int(8, 5), 0.1, 0.4, i))
})

simulation_results_gooogleplus <- measures.summary(goooglePlus, 200, data.list, "grLasso", 5, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(simulation_results_gooogleplus)
simulation_results_gooogle <- measures.summary(gooogle, 200, data.list, "grLasso", 5, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)))
print(simulation_results_gooogle)
