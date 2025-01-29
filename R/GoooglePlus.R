library(pscl)
library(vscDebugger)
library(glmnet)
library(mpath)
library(rockchalk)
library(Gooogle)
library(outliers)
library(forecast)
library(gamlss.dist)

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
    evcov <- eigen(vcov, symmetric = TRUE)
    # vcov is symmetric square matrix therefore the matrix of its eigenvectors should be orthogonal ie AT = A^-1
    inv_sqrt_diag <- diag(vapply(evcov$values, function(num) {
        num ^ (-1 / 2)
    }, numeric(1)), nrow = nrow(evcov$vectors), ncol = ncol(evcov$vectors))
    pseudo_factor <- evcov$vectors %*% inv_sqrt_diag %*% t(evcov$vectors)

    pseudo_X <- pseudo_factor
    pseudo_Y <- as.vector(pseudo_factor %*% b2.mle)

    # standardise pseudo_X (not pseudo_Y)
    # scale only and do not centre as we do not want to produce another intercept
    X_std <- standardise(pseudo_X)
    pseudo_X <- X_std$scaled_X
    X_scale <- X_std$s
    fit.coefficients <- fit.coefficients / sqrt(mean(fit.coefficients ^ 2))

    # reorder the group index and the covariates
    group <- c(0, group.x, 0, group.z)

    # initialise pseudo_X reordered by group
    group_ordered_indices <- order(group)
    # get the indices for putting into unordered form again (original form)
    unordering_indices <- numeric(length(group))
    for (i in 1:length(group_ordered_indices)) {
        unordering_indices[group_ordered_indices[i]] <- i
    }
    pseudo_X <- pseudo_X[, group_ordered_indices]
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

    # Step 1: De-orthonormalise and unorder coefficients and de-scale
    coefficients_list <- lapply(orthonormalised_coefficients_list, function(coeff) {
        return((deorthonormalise_coefficients(
            coeff,
            orthonormalisation_factors,
            group_start_indices
        )[unordering_indices]) / X_scale)
    })

    mat_X <- as.matrix(cbind(1, X))
    mat_Z <- as.matrix(cbind(1, Z))
    vect_y <- as.numeric(y[, 1])

    # Step 2: Find bic and optimise
    #coefficient_log_likelihood <- vapply(coefficients_list, function(coeff) {
    #    return(
    #        zip_log_likelihood(
    #            mat_X,
    #            mat_Z,
    #            coeff[1:(p + 1)],
    #            coeff[(p + 2): (p + q + 2)],
    #            mat_y
    #        )
    #    )
    #}, numeric(1))

    coefficient_log_likelihood <- vapply(coefficients_list, function(coeff) {
        return(ll.func(coeff[1:(p + 1)], coeff[(p + 2): (p + q + 2)], vect_y, mat_X, mat_Z))
    }, numeric(1))

    # get number of model parameters (non-zero coefficients) and exclude intercepts - we have two intercepts
    # at the moment, this may be overstated => Zeng's approach of e^lambda may be better.
    df <- vapply(coefficients_list, function(x) { return(sum(x != 0) - 2) }, numeric(1))
    coefficients_bic <- vapply(1:length(coefficient_log_likelihood), function(i) {
        return(df[i] * log(nrow(X)) - 2 * coefficient_log_likelihood[i])
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
                log_likelihood = coefficient_log_likelihood[min_bic_index],
                bic = optimal_bic,
                coefficients = opt_params,
                opt_lambda = optimal_lambda,
                params = coefficients_list))
}

orthonormalise_full_matrix <- function(all_groups_matrix, group_start_indices) {
    n <- nrow(all_groups_matrix)
    orthonormalised_matrix <- matrix(0, nrow = n, ncol = ncol(all_groups_matrix))
    matrix_factors <- vector("list", length(group_start_indices))
    intercept_factors <- numeric(group_start_indices[2] - group_start_indices[1])
    for (i in group_start_indices[1]:(group_start_indices[2] - 1)) {
        orth_intercept <- svd_intercept_vector(all_groups_matrix[, i], n)
        orthonormalised_matrix[, i] <- orth_intercept[[1]]
        intercept_factors[i] <- orth_intercept[[2]]
    }
    matrix_factors[[1]] <- intercept_factors
    for (i in 2:(length(group_start_indices) - 1)) {
        # get the matrix:
        group_matrix <- all_groups_matrix[, group_start_indices[i]:(group_start_indices[i + 1] - 1)]
        orth_group_matrix_and_factor <- svd_group_matrix(group_matrix, n)
        orthonormalised_group_matrix <- orth_group_matrix_and_factor[[1]]
        matrix_factors[[i]] <- orth_group_matrix_and_factor[[2]]
        orthonormalised_matrix[, group_start_indices[i]:(group_start_indices[i + 1] - 1)] <- orthonormalised_group_matrix
    }
    return(list(orthonormalised_matrix, matrix_factors))
}

svd_group_matrix <- function(matrix, n) {
    svdresult <- svd(matrix, nu = 0)
    # r <- which(svdresult$d > 1e-10)
    matrix_factor <- sqrt(n) * svdresult$v %*% diag(vapply(svdresult$d, function(num) {
        return(1 / num)
    }, numeric(1)))
    matrix_and_factor <- vector("list", 2)
    matrix_and_factor[[1]] <- matrix %*% matrix_factor
    matrix_and_factor[[2]] <- matrix_factor
    return(matrix_and_factor)
}

svd_intercept_vector <- function(vect, n) {
    svdresult <- svd(vect, nu = 0)
    vfactor <- sqrt(n) * svdresult$v * 1 / svdresult$d
    return(list(
        vect %*% vfactor,
        vfactor
    ))
}

orthogonalize <- function(X, group) {
    n <- nrow(X)
    J <- max(group)
    T <- vector("list", J)
    XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
    XX[, which(group==0)] <- X[, which(group==0)]
    for (j in seq_along(integer(J))) {
        ind <- which(group==j)
        if (length(ind)==0) next
        SVD <- svd(X[, ind, drop=FALSE], nu=0)
        r <- which(SVD$d > 1e-10)
        T[[j]] <- sweep(SVD$v[, r, drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
        XX[, ind[r]] <- X[, ind] %*% T[[j]]
    }
    # do this after function calls to orthonormalise group matrices
    nz <- !apply(XX==0, 2, all)
    XX <- XX[, nz, drop=FALSE]
    attr(XX, "T") <- T
    attr(XX, "group") <- group[nz]
    XX
}

unorthogonalize <- function(b, XX, group, intercept=TRUE) {
    require(Matrix)
    ind <- !sapply(attr(XX, "T"), is.null)
    T <- bdiag(attr(XX, "T")[ind])
    val <- as.matrix(T %*% b)
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

orthonormalise_intercept_vector <- function(vect) {
    factor <- sqrt(length(vect) / sum(vect ^ 2))
    return(list(
        factor * vect,
        factor
    ))
}

deorthonormalise_coefficients <- function(coeff, matrix_factors, group_start_indices) {
    deorth_coeff <- numeric(length(coeff))
    for (i in group_start_indices[1]:(group_start_indices[2] - 1)) {
        deorth_coeff[i] <- matrix_factors[[1]][i] * coeff[i]
    }
    for (i in 2:(length(group_start_indices) - 1)) {
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
        standard_deviation_y <- sqrt(sum(y ^ 2) / nrow(group_mat_orth))
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
                X <- group_mat_orth[, group_start_indices[j]] # will not be changed
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
                X <- group_mat_orth[, group_start_indices[j]:(group_start_indices[j + 1] - 1)] # will not be changed
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

standardise <- function(X) {
    # There is no consideration of grouping when performing feature-level standardisation
    # pseudo_X is already a matrix
    p <- ncol(X)
    s <- numeric(p)

    for (j in seq_len(p)) {
        # Do not centre

        # Scale
        s[j] <- sqrt(mean(X[, j]^2))
        X[, j] <- X[, j] / s[j]
    }

    # Return list
    res <- list(scaled_X = X, s = s)
    return(res)
}

unstandardise_coefficients <- function(coeff, stddev) {
    beta <- numeric(length(coeff))
    beta <- coeff / stddev
    return(beta)
}

ll.func <- function(beta.count, beta.zero, y, X, Z, dist)
{
    if (is.null(Z))
    {
        zgam <- rep(beta.zero, length(y))
    } else {
        zgam <- Z %*% beta.zero
        zgam <- as.numeric(zgam[, 1])
    }
    pzero <- exp(zgam) / (1 + exp(zgam))

    xbet <- X %*% beta.count
    xbet <- as.numeric(xbet[, 1])
    mu <- exp(xbet)

    ll <- try(sum(dZIP(y, mu = mu, sigma = pzero, log = TRUE)), silent = TRUE)

    # error-handling
    if (class(ll) == "try-error")
    {
        ll <- NA
    }
    return(ll)
}

# Improvements:

# Standardise the whole vector before using group standardisation:
# Find a way to do this without requiring an intercept
# We have centering and scaling
# What if we only scale and do not centre

# Try to remove intercepts as they do
output <- gen_zip_data(200, 50, rep.int(8, 5), 0.1, 0.4, 200)
data <- output$data
yvar <- output$yvar
xvars <- output$xvars
zvars <- output$zvars
sim_result_grLasso <- goooglePlus(data, xvars, zvars, yvar, c(rep(1, 8), rep(2, 8), rep(3, 8), rep(4, 8), rep(5, 8)), penalty = "grLasso", lambda_min = 0)
print(sim_result_grLasso$coefficients)
print(sim_result_grLasso$params[[100]])
