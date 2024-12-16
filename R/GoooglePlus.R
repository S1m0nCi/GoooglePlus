install.packages("zeroinfl")
library(zeroinfl)

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
    lambda.min = ifelse((nrow(data[, unique(c(xvars, zvars))]) > ncol(data[, unique(c(xvars, zvars))])), 1e-4, .05),
    lambda.max,
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
    n <- nrow(data)

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
    # the reason for doing all of the above, as we see in the 'Insurance' example of calling gooogle, is because we usually will set zvars = xvars, but in general the zeroinf
    # explanatory variables and the count explanatory variables may not be in the same groups.
    # get p and q for use later with the intercepts.
    # vcov matrix also includes for the intercepts
    p <- length(xvars)
    q <- length(zvars)
    # use zeroinfl to fit a model, which is part of pscl: this requires "|" in the formula as has been set above - that is the reason why we did it all
    fit.zero <- zeroinfl(fit.formula, dist = dist, data = data)
    b2.mle <- c(fit.zero$coefficients$count[-1], fit.zero$coefficients$zero[-1])
    # INVESTIGATE: create vector of coefficients with the intercept coefficients set to zero
    fit.coefficients <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)
    zeroinf_residuals <- fit.zero$residuals
    # get the covariance matrix needed for transforming the data
    vcov <- fit.zero$vcov #  TODO: this needs to have the correct parts removed, ie the intercepts need to be removed.
    # currently vcov is a square p + q + 2 matrix.
    # how do we make the intercepts zero and allow them to be updated? The gcd is only done for groups and we are not counting the intercepts as groups.

    p <- length(xvars)

    # possible gimmicks start here up to transformed y and x
    # we will need to do the pseudo data as done by Zeng

    # get the pseudo data factor sigma to the -1/2
    evcov <- eigen(vcov)
    # vcov is a square, symmetric matrix
    if (det(evcov) > 0) {  # in case vcov is not pd add small values to the diagonal ie makePD
        pseudo_factor <- evcov$vectors %*% diag(1 / sqrt(abs(evcov$values))) %*% t(evcov$vectors)
    } else {
        pseudo_factor <- evcov$vectors %*% diag(1 / sqrt(makePD(evcov))) %*% t(evcov$vectors)
    }

    pseudo_X <- pseudo_factor
    pseudo_Y <- pseudo_factor %*% b2.mle

    # 2 by 2 matrix
    # sigma.11 <- vcov[c(1, p + 2), c(1, p + 2)]
    # 2 by p matrix
    # sigma.12 <- vcov[c(1, p + 2), -c(1, p + 2)]
    # p by p matrix
    # sigma.22 <- vcov[-c(1, p + 2), -c(1, p + 2)]

    # vcov.bar <- sigma.22 - t(sigma.12) %*% ginv(sigma.11) %*% sigma.12
    # e <- eigen(vcov.bar)

    # if (det(vcov.bar) > 0) {  # in case vcov is not pd add small values to the diagonal ie makePD
    #      cov.star <- e$vectors %*% diag(1 / sqrt(abs(e$values))) %*% t(e$vectors)
    # } else {
    #     cov.star <- e$vectors %*% diag(1 / sqrt(makePD(vcov.bar))) %*% t(e$vectors)
    # }

    # y.star <- cov.star %*% b2.mle  # transformed y
    # cov.star <- data.frame(cov.star)  # scaled x matrix
    # names(cov.star) <- c(xvars, zvars)

    # reorder the group index and the covariates
    # TODO: Change to put all intercepts in one group numbered 0.
    group <- c(group.x, 0, group.z, 0)

    unique_groups <- unique(group)
    # initialise pseudo_X reordered by group
    group_ordered_indices <- order(group)
    # get the indices for putting into unordered form again (original form)
    unordering_indices <- numeric(length(group))
    for (i in 1:length(group_ordered_indices)) {
        unordering_indices[group_ordered_indices[i]] <- i
    }
    pseudo_X <- pseudo_X[, group_ordered_indices]
    fit.coefficients <- fit.coefficients[group_ordered_indices]
    zeroinf_residuals <- zeroinf_residuals[group_ordered_indices]

    # pseudo_X.reordered <- rep(0, dim(cov.star)[1])
    # for (i in 1:length(uniqu_groups)) {
    #     indices <- which(group == unique_groups[i])
    #     cov.star.reordered <- cbind(cov.star.reordered, cov.star[, indices])
    # }
    # cov.star.reordered <- cov.star.reordered[, -1]
    group <- sort(group)

    # scaling and group sorting has now already been done
    # now we implement the group coordinate descent algorithm
    # according to Zeng, we do not need to transform back to non-LSA, so non-pseudo data
    # however we will need to orthogonalise and then transform back

    # get the group matrices
    # use 'match'
    group_start_indices <- c(match(unique_groups, group), length(group) + 1)
    # group_matrices <- vector("list", length(unique_groups))
    penalty_factors <- numeric(length(unique_groups))
    # the intercepts are to be unpenalised - setting penalisation to 0.
    penalty_factors[0] <- 0
    group_initial_coefficients <- vector("list", length(unique_groups))

    # we want penalty factors of zero for the intercepts, so the last two elements of the group
    for (i in 2:length(group_start_indices) - 1) {
        # group_matrices[[i]] <- pseudo_X[, group_start_indices[i]:group_start_indices[i + 1] - 1]
        penalty_factors[i] <- sqrt(group_start_indices[i + 1] - group_start_indices[i])
        group_initial_coefficients[[i]] <- fit.coefficients[group_start_indices[i]:group_start_indices[i + 1] - 1]
    }
    # group_matrices[length(unique_groups)] <- pseudo_X[, group_start_indices[length(unique_groups)]:ncol(pseudo_X)]

    group_initial_coefficients[length(unique_groups)] <- fit.coefficients[group_start_indices[length(unique_groups)]:length(fit.coefficients)]

    # orthonormalise the group matrices
    group_matrices_orth_full <- orthonormalise_full_matrix(pseudo_X)
    group_matrix_orthonormal <- group_matrices_orth_full[1]
    orthonormalisation_factors <- group_matrices_orth_full[2]

    # Cross Validation - one parameter only, which is lambda
    if (missing(lambda)) {
        if (missing(lambda.max)) {
            intercepts_fit <- zeroinfl(y ~ 1 | 1) # will give matrix of 2 cols, instead of vector
            zmax <- maxgrad(pseudo_X, intercepts_fit$residuals, group_start_indices, penalty_factors) / nrow(pseudo_X)
            lambda.max <- zmax
        }
        if (log.lambda) {
            # creating in descending order
            if (lambda.min == 0) {
                lambda <- c(exp(seq(log(lambda.max), log(0.001 * lambda.max), length = nlambda - 1)), 0)
            } else {
                lambda <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), length = nlambda))
            }
        }
        if (lambda.min == 0) {
            lambda <- c(seq(lambda.max, 0.001 * lambda.max, length = nlambda - 1), 0)
        } else {
            lambda <- seq(lambda.max, lambda.min * lambda.max, length = nlambda)
        }
    }

    orthonormalised_coefficients_list <- group_coordinate_descent(
        group_matrix_orthonormal,
        lambda,
        penalty_factors,
        group_initial_coefficients,
        zeroinf_residuals,
        eps,
        group_start_indices,
        pseudo_Y
    )
    # TODO: Post-processing
    # Step 1: De-orthonormalise and unorder coefficients
    deorthonormalised_coefficients_list <- vector("list", length(orthonormalised_coefficients_list))
    for (i in 1:orthonormalised_coefficients_list) {
        deorthonormalised_coefficients_list[[i]] <- deorthonormalise_coefficients(orthonormalised_coefficients_list[i], orthonormalisation_factors, group_start_indices)[unordering_indices]
    }

    # Step 2: de-scale - this may not be needed as we used vcov to generate pseudo data, as according to Zeng.

    # Step 3: Find bic and optimise
    loss_of_coefficients <- lapply(deorthonormalised_coefficients_list, function(coefficients) {
        return(mean_squared_error(as.matrix(X), coefficients, as.vector(y)))
    })

    # get number of model parameters (non-zero coefficients) and exclude intercept
    df <- lapply(deorthonormalised_coefficients_list, 2, function(x) { return(sum(x != 0)) }) - 1
    coefficients_bic <- numeric(deorthonormalised_coefficients_list)
    for (i in 1:length(deorthonormalised_coefficients_list)) {
        coefficients_bic[i] <- df[[i]] * log(nrow(X)) + nrow(X) * log(loss_of_coefficients[[i]] / nrow(X))
    }

    min_bic_index <- which.min(coefficients_bic)
    optimal_coefficients <- deorthonormalised_coefficients_list[[min_bic_index]]
    optimal_lambda <- lambda[i]
    optimal_bic <- coefficients_bic[min_bic_index]

    # TO DO: consider un-ordering the coefficients: after de-orthonormalising as we used the ordered/grouped ones there
    return(list(
                params = deorthonormalised_coefficients_list,
                group = group,
                lambda = lambda,
                df = df,
                loss = loss_of_coefficients,
                bic = coefficients_bic,
                penalty = penalty,
                n = nrow(X),
                # iter = iterations,
                opt_params = optimal_coefficients,
                opt_lambda = optimal_lambda,
                opt_bic = optimal_bic))
}

orthonormalise_full_matrix <- function(all_groups_matrix, group_start_indices) {
    num_cols <- ncol(all_groups_matrix)
    orthonormalised_matrix <- matrix(, nrow = nrow(all_groups_matrix), ncol = num_cols)
    matrix_factors <- vector("list", length(group_start_indices))
    for (i in 1:length(group_start_indices) - 1) {
        # get the matrix:
        group_matrix <- all_groups_matrix[, group_start_indices[i]:group_start_indices[i + 1] - 1]
        orth_group_matrix_and_factor <- orthonormalise_group_matrix(group_matrix)
        orthonormalised_group_matrix <- orth_group_matrix_and_factor[1]
        matrix_factors[[i]] <- orth_group_matrix_and_factor[2]
        orthonormalised_matrix[, group_start_indices[i]:group_start_indices[i + 1] - 1] <- orthonormalised_group_matrix
    }
    return(c(orthonormalised_group_matrix, matrix_factors))
}

# vapply is faster then lapply

orthonormalise_group_matrix <- function(matrix) {
    gram <- (1 / nrow(matrix)) * t(matrix) %*% matrix
    eigenresult <- eigen(gram)
    eigenvector_matrix <- eigenresult$vectors
    inv_diag_matrix <- diag(sapply(eigenresult$values), function(num) {
        num ^ (-1 / 2)
    })
    matrix_factor <- eigenvector_matrix * inv_diag_matrix
    # keep the matrix factor for transformation back later on.
    return(c(matrix %*% matrix_factor, matrix_factor))
}

deorthonormalise_coefficients <- function(coeff, matrix_factors, group_start_indices) {
    deorth_coeff <- numeric(length(coeff))
    for (i in 1:length(group_start_indices - 1)) {
        deorth_coeff[[group_start_indices[i]:group_start_indices[i + 1] - 1]] <- matrix_factors[i] %*% coeff[[group_start_indices[i]:group_start_indices[i + 1] - 1]]
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
    gamma = 0,
    penalty = "lasso"
) {
    # we need to output all of the coefficients from the different lambdae
    # we then calculate loss etc using bic later
    # also need to worry about how to reorder
    # 'coefficients' below may need to be a vector instead of a list
    coefficients_list <- vector("list", length(lambda))

    for (i in 1:length(lambda)) {
        lam <- lambda[i]

        # do not keep coefficients in group structure
        coefficients <- numeric(ncol(group_mat_orth))
        standard_deviation_y <- sqrt(sum(y ** 2) / nrow(group_mat_orth))
        convergence_threshold <- eps * standard_deviation_y
        penalty_factors <- penalty_factors * lam
        # do the intercepts first. We treat them as separate, not as one group. We just 'store' them as one group.
        for (j in group_start_indices[1]:group_start_indices[2] - 1) {
            iteration <- 1
            lambda <- penalty_factors[1]
            b <- initial_b[1][j]
            r <- initial_residuals
            X <- group_mat_orth[, group_start_indices[j]] # will not be changed
            while (iteration < max_iterations) {
                z <- t(X) %*% r + b
                new_b <- vector_soft_threshold(z, lambda_i)
                b_change <- new_b - b
                r <- r - X %*% (b_change) # is there an error in the paper by Breheny, which says to use t(X)?
                max_change <- max(b_change)
                if (max_change < convergence_threshold) {
                    break
                }
                iteration <- iteration + 1
            }
            coefficients[j] <- b
        }

        for (j in 2:length(group_start_indices) - 1) {
            iteration <- 1
            lambda_i <- penalty_factors[j]
            b <- initial_b[j]
            r <- initial_residuals
            X <- group_mat_orth[, group_start_indices[j]:group_start_indices[j + 1] - 1] # will not be changed
            while (iteration < max_iterations) {
                z <- t(X) %*% r + b
                new_b <- vector_soft_threshold(z, lambda_i)
                b_change <- new_b - b
                r <- r - X %*% (b_change) # is there an error in the paper by Breheny, which says to use t(X)?
                max_change <- max(b_change)
                if (max_change < convergence_threshold) {
                    break
                }
                iteration <- iteration + 1
            }
            coefficients[group_start_indices[j]:group_start_indices[j + 1] - 1] <- b
        }
        coefficients_list[i] <- coefficients
    }
    return(coefficients)
}
# TO DO: Integrate the different penalty functions with the group coordinate descent.
# would be useful if they had types
# lasso for now, group MCP and SCAD for later
lasso_penalty <- function(vect) {
    return(sqrt(sum(vect^2)))
}

mcp_penalty <- function() {}
scad_penalty <- function() {}

soft_threshold <- function(num, lambda) {
    if (num > lambda) {
        return(num - lambda)
    } else if (num < - lambda) {
        return(lambda - num)
    }
    return(0)
}

vector_soft_threshold <- function(vect, lambda) {
    norm <- sqrt(sum(vect^2))
    return(soft_threshold(norm, lambda) * vect / norm)
}

maxgrad <- function(pseudo_X, residuals, group_start_indices, penalty_factors) {
    # XX_orth, r, K1, group_multipliers
    zmax <- 0

    # ignore the intercepts here as they are not affected by lambda, as they are not penalised
    for (i in 2:length(group_start_indices) - 1) {
        # get group size
        group_length <- group_start_indices[i + 1] - group_start_indices[i]
        Z <- numeric(group_length)

        for (j in (group_start_indices[i] + 1):(group_start_indices[i + 1])) {
            # dot product of column and residuals
            Z[j - group_start_indices[i]] <- sum(pseudo_X[, j] * residuals)
        }
        z <- sqrt(sum(Z^2)) / penalty_factors[i]
        zmax <- max(zmax, z)
    }
    return(zmax)
}

mean_squared_error <- function(pseudo_X, coeff, actual_obs) {
    return(sum((pseudo_X %*% coeff - actual_obs)^2))
}
