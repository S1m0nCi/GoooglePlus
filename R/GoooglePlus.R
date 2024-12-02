install.packages("zeroinfl")

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

    # use zeroinfl to fit a model, which is part of pscl: this requires "|" in the formula as has been set above - that is the reason why we did it all
    fit.zero <- zeroinfl(fit.formula, dist = dist, data = data)
    b2.mle <- c(fit.zero$coefficients$count[-1], fit.zero$coefficients$zero[-1])
    fit.coefficients <- fit.zero$coefficients

    # get the covariance matrix needed for transforming the data
    vcov <- fit.zero$vcov
    p <- length(xvars)

    # possible gimmicks start here up to transformed y and x
    # we will need to do the pseudo data as done by Zeng

    # get the pseudo data factor sigma to the -1/2
    evcov <- eigen(vcov)
    # vcov is a square, symmetric matrix
    if (det(evcov) > 0) {  # in case vcov is not pd add small values to the diagonal ie makePD
        pseudo_factor <- evcov$vectors %*% diag(1 / sqrt(abs(evcov$values))) %*% t(evcov$vectors)
    } else {
        pseudo_factor <- evcov$vectors %*% diag(1 / sqrt(makePD(evcov))) %*% t(e$vectors)
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
    group <- c(group.x, group.z)
    unique_groups <- unique(group)
    # initialise pseudo_X reordered by group
    group_ordered_indices <- order(group)
    pseudo_X <- pseudo_X[, group_ordered_indices]

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
    group_start_indices <- match(unique_groups, group)
    group_matrices <- numeric(length(unique_groups))
    penalty_factors <- numeric(length(unique_groups))
    group_initial_coefficients <- numeric(length(unique_groups))
    zeroinf_residuals <- fit.zero$residuals

    for (i in 1:length(group_start_indices)-1) {
        group_matrices[i] <- pseudo_X[, group_start_indices[i]:group_start_indices[i + 1] - 1]
        penalty_factors[i] <- sqrt(group_start_indices[i + 1] - group_start_indices[i])
        group_initial_coefficients[i] <- fit.coefficients[group_start_indices[i]:group_start_indices[i + 1] - 1]
    }
    group_matrices[length(unique_groups)] <- pseudo_X[, group_start_indices[length(unique_groups)]:ncol(pseudo_X)]
    penalty_factors[length(unique_groups)] <- sqrt(ncol(pseudo_X) - group_start_indices[length(unique_groups) + 1])
    group_initial_coefficients[length(unique_groups)] <- fit.coefficients[group_start_indices[length(unique_groups)]:length(fit.coefficients)]

    # orthonormalise the group matrices
    group_matrices_orth_full <- sapply(group_matrices, orthonormalise)
    group_matrices_orthonormal <- sapply(group_matrices_orth_full, function(pair) {
        return(pair[2])
    })
    orthonormalisation_factors <- sapply(group_matrices_orth_full, function(pair) {
        return(pair[1])
    })

    orthonormalised_coefficients <- group_coordinate_descent(
        group_matrices_orthonormal,
        lambda,
        penalty_factors,
        group_initial_coefficients,
        zeroinf_residuals,
        eps,
        pseudo_Y
    )
    
    final_coefficients <- numeric(length(orthonormalised_coefficients))
    for (i in 1:length(orthonormalised_coefficients)) {
        final_coefficients[i] <- deorthonormalise_coeffiecients(orthonormalised_coefficients[i], orthonormalisation_factors[i])
    }

    return(final_coefficients)



    # TO DO: Add the intercept update in group_coordinate_descent and return the intercepts in the function
    # do not penalise the intercepts

}

orthonormalise <- function(matrix) {
    gram <- (1 / nrow(matrix)) * t(matrix) %*% matrix
    eigenresult <- eigen(gram)
    eigenvector_matrix <- eigenresult$vectors
    inv_diag_matrix <- diag(sapply(eigenresult$values), function(num) {
        num ^ (-1 / 2)
    }) # vapply faster
    matrix_factor <- eigenvector_matrix * inv_diag_matrix
    # keep the matrix factor for transformation back later on.
    return(c(matrix %*% matrix_factor, matrix_factor))
}

deorthonormalise_coeffiecients <- function(coeff, matrix_factor) {
    return(matrix_factor %*% coeff)
}

l2_distance <- function(u, v) {
    return(sqrt(sum((u - v) ^ 2)))
}

# TO DO: Implement group coordinate descent
group_coordinate_descent <- function(
    group_matrices, # the list of matrices for groups
    lambda,
    penalty_factors,
    initial_b,
    initial_residuals,
    eps,
    y,
    max_iterations = 1000,
    gamma = 0,
    penalty = "lasso"
) {
    coefficients <- numeric(length(group_matrices))
    standard_deviation_y <- sqrt(sum(y ** 2) / sum(sapply(group_matrices, function(mat) {return(nrow(mat))})))
    convergence_threshold <- eps * standard_deviation_y
    penalty_factors <- penalty_factors * lambda

    for (i in 1:length(group_matrices)) {
        iteration <- 1
        lambda_i <- penalty_factors[i]
        b <- initial_b[i]
        r <- initial_residuals[i]
        X <- group_matrices[i] # will not be changed
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
        coefficients[i] <- b
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
