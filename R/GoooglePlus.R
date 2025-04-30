library(pscl)
library(vscDebugger)
library(gamlss.dist)
library(roxygen2)
library(stats)
library(numDeriv)

library(rockchalk)
library(mpath)

dyn.load("C:\\Users\\gbsim\\Documents\\R-scripts\\GoooglePlus\\src\\gcd.dll")

#' A group regularized fit to the zero inflated count data.
#'
#' @description Fit zero inflated count data with a group regularization algorithm.
#'
#' @usage gooogleplus(data,xvars,zvars,yvar,group=1:ncol(data),samegrp.overlap=T,penalty=c("grLasso", "grMCP", "grSCAD", "gBridge"),dist=c("poisson","negbin"), nlambda=100, lambda,lambda_min=ifelse((nrow(data[,unique(c(xvars,zvars))])>ncol(data[,unique(c(xvars,zvars))])),1e-4,.05),lambda_max, crit="BIC",alpha=1, eps=.001, max_iter=1000, gamma=ifelse(penalty=="gBridge",0.5,ifelse(penalty == "grSCAD", 4, 3)))
#'
#'
#' @param data The data frame or matrix consisting of outcome and predictors.
#' @param xvars The vector of variable names to be included in count model.
#' @param zvars The vector of variable names for excess zero model.
#' @param yvar The outcome variable name.
#' @param group The vector of integers describing the grouping of the coefficients. For greatest efficiency and least ambiguity, it is best if group is a vector of consecutive integers. If there are coefficientss to be included in the model without being penalized, assign them to group 0 (or "0").
#' @param samegrp.overlap A logical argument. If TRUE (default) same grouping indices will be assigned to shared predictors in the count and degenerate distribution.
#' @param penalty The penalty to be applied in the model. For group level selection, one of "grLasso", "grMCP" or "grSCAD". For bi-level selection "gBridge" can be specified.
#' @param dist The distribution for count model - "poisson" for poisson or "negbin" for negative binomial.
#' @param nlambda The number of lambda values. Default is 100.
#' @param lambda A user specified sequence of lambda values.
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .0001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param lambda_max The maximum value for lambda (only needed for gBridge penalty).
#' @param crit The selection criteria for the best model. It can either be "AIC" or \code{BIC} (default).
#' @param alpha The tuning parameter for the balance between the group penalty and the L2 penalty, as in grpreg. Default value is 1.
#' @param eps The convergence threshhold.
#' @param max_iter Maximum number of iterations allowed.
#' @param gamma Tuning parameter of group MCP/SCAD. Default is 3 for MCP and 4 for SCAD.
#'
#' @details The algorithm fits zero inflated count data to conduct variable selection in the presence of intrinsic grouping structure in the predictor set. Group wise penalties are considered for both count and zero abundance part of the mixture model where the likelihood is optimized using group level co-ordinate descent algorithms.
#'
#' @return A list containing the following components is returned
#' \item{log_likelihood}{The log-likelihood of the selected model.}
#' \item{bic}{The BIC of the selected model.}
#' \item{coefficients}{A list with two sets of coefficients corresponding to count and zero inflation parts of the mixture model.}
#' \item{opt_lambda}{The lambda value chosen for the model.}
#'
#' @export
#' @examples
#' \dontrun{
#' ## Auto Insurance Claim Data
#' library(HDtweedie)
#' data("auto")
#' y<-auto$y
#' y<-round(y)
#' x<-auto$x
#' data<-cbind.data.frame(y,x)
#' group=c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)
#' yvar<-names(data)[1]
#' xvars<-names(data)[-1]
#' zvars<-xvars
#'
#' ## ZIP regression
#' fit.poisson<-gooogleplus(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=TRUE,dist="poisson",penalty="grLASSO")
#' fit.poisson$bic
#'
#' ## ZINB regression
#' fit.negbin<-gooogleplus(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=TRUE,dist="negbin",penalty="grLASSO")
#' fit.negbin$bic
#' }
#'
#' @importFrom pscl zeroinfl
#' @importFrom stats glm
#' @importFrom gamlss.dist dZIP dZINBI
#'
gooogleplus <-  function(
    data,
    xvars,
    zvars,
    yvar,
    group = seq_len(ncol(data)),
    samegrp.overlap = TRUE,
    penalty = c("grLasso", "grMCP", "grSCAD", "grALasso"),
    dist = c("poisson", "negbin"),
    nlambda = 100,
    lambda,
    lambda_min = ifelse((nrow(data[, unique(c(xvars, zvars))]) > ncol(data[, unique(c(xvars, zvars))])), 1e-4, .05),
    lambda_max,
    log.lambda = TRUE,
    crit = "BIC",
    alpha = 1,
    eps = .001,
    max_iter = 1000,
    gamma = ifelse(penalty == "grSCAD", 4, 3),
    runtime = c("C", "R", "compare")
)
{
    runtime <- match.arg(runtime)
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

        # the below code just ensures no overlap of X and Z group numbers: they are treated as completely separate groups:
        # hence we add the max of groups x to ensure none are shared
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

    a <- 1
    if (dist == "negbin") {
        a <- 1 / fit.zero$theta
    }

    # keep the intercepts in the mle and the final data, which we want to have p + q + 2 columns [Zeng]
    # different to Gooogle where the approximation is used.
    b2.mle <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)

    fit.coefficients <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)
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
    fit.coefficients <- fit.coefficients * X_scale

    mat_X <- as.matrix(cbind(1, X))
    mat_Z <- as.matrix(cbind(1, Z))
    vect_y <- as.numeric(y[, 1])

    if (penalty == "BAR") {
        if (missing(lambda)) {
            if (missing(lambda_max)) {
                colnames(pseudo_X) <- names(fit.coefficients)
                bar_intercepts_fit <- glm(
                    as.formula("y ~ 0 + count.intercept + zero.intercept"),
                    data = cbind(data.frame(y = unlist(pseudo_Y)), pseudo_X)
                )
                # now deorthonormalise those
                bar_intercepts_coeff <- rep(0, p + q + 2)
                bar_intercepts_coeff[1] <- bar_intercepts_fit$coefficients["count.intercept"] / X_scale[1]
                bar_intercepts_coeff[p + 2] <- bar_intercepts_fit$coefficients["zero.intercept"] /  X_scale[p + 2]
                lambda_max <- maxgrad_bar(mat_X, mat_Z, vect_y, bar_intercepts_coeff, p, q)
            }
            lambda <- create_lambda(lambda_max, lambda_min, log.lambda, nlambda)
        }
        else {
            nlambda <- length(lambda)
        }
        coefficients_list <- cycBAR(mat_X, mat_Z, vect_y, pseudo_X, pseudo_Y, X_scale, lambda, p, q, max_iter, eps)
        print(coefficients_list)
    } else {
        # reorder the group index and the covariates
        group <- c(0, group.x, 0, group.z)

        # initialise pseudo_X reordered by group
        group_ordered_indices <- order(group)
        # get the indices for putting into unordered form again (original form)
        unordering_indices <- numeric(length(group))
        for (i in seq_along(group_ordered_indices)) {
            unordering_indices[group_ordered_indices[i]] <- i
        }
        pseudo_X <- pseudo_X[, group_ordered_indices]
        fit.coefficients <- fit.coefficients[group_ordered_indices]

        group <- group[group_ordered_indices]
        unique_groups <- unique(group)

        # get the group matrices
        # use 'match' and convert to integer for C conversion
        group_start_indices <- as.integer(c(match(unique_groups, group), length(group) + 1))
        penalty_factors <- numeric(length(unique_groups))
        # the intercepts are to be unpenalised - setting penalisation to 0.
        penalty_factors[1] <- 0

        # we want penalty factors of zero for the intercepts, so the last two elements of the group
        for (i in 2:(length(group_start_indices) - 1)) {
            # group_matrices[[i]] <- pseudo_X[, group_start_indices[i]:group_start_indices[i + 1] - 1]
            penalty_factors[i] <- sqrt(group_start_indices[i + 1] - group_start_indices[i])
        }

        if (penalty == "grALASSO") {
            # length of penalty factors and length of group norms should be the same
            for (i in 2:(length(group_start_indices) - 1)) {
                group_norm <- sqrt(sum(fit.coefficients[group_start_indices[i]:group_start_indices[i + 1]]))
                if (group_norm == 0) {
                    # avoid zero division
                    group_norm <- .Machine$double.eps
                }
                penalty_factors[i] <- penalty_factors[i] / group_norm
            }
        }

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
                    as.formula(paste("y ~ 0 + count.intercept + zero.intercept")),
                    data = cbind(data.frame(y = unlist(pseudo_Y)), group_matrix_df)
                )
                # will give matrix of 2 cols, instead of vector
                # make a matrix of X, Z and rearrange for groups as necessary
                zmax <- maxgrad(group_matrix_orthonormal, intercepts_fit$residuals, group_start_indices, penalty_factors) / nrow(group_matrix_orthonormal)
                lambda_max <- zmax
            }
            lambda <- create_lambda(lambda_max, lambda_min, log.lambda, nlambda)
        } else {
            nlambda <- length(lambda)
        }

        if (runtime == "R" || runtime == "compare") {

            orthonormalised_coefficients_list <- group_coordinate_descent(
                group_matrix_orthonormal,
                lambda,
                penalty_factors,
                fit.coefficients,
                eps,
                group_start_indices,
                pseudo_Y,
                penalty = penalty,
                max_iterations = max_iter
            )

            coefficients_list <- lapply(orthonormalised_coefficients_list, function(coeff) {
                return(deorthonormalise_coefficients(
                    coeff,
                    orthonormalisation_factors,
                    group_start_indices
                )[unordering_indices] / X_scale)
            })
        }
        if (runtime == "C" || runtime == "compare") {
            orthonormalised_coefficients_cm <- .Call(
                "groupCoordinateDescent",
                group_matrix_orthonormal,
                lambda,
                penalty_factors,
                fit.coefficients,
                eps,
                group_start_indices,
                pseudo_Y,
                gamma,
                penalty,
                as.integer(max_iter)
            )

            group_start_indices <- group_start_indices + rep(1, length(group_start_indices))

            orthonormalised_coefficients <- matrix(orthonormalised_coefficients_cm, nrow = p + q + 2, ncol = nlambda)
            # Step 1: De-orthonormalise and unorder coefficients
            coefficients_list <- apply(orthonormalised_coefficients, 2, function(coeff) {
                return(deorthonormalise_coefficients(
                    coeff,
                    orthonormalisation_factors,
                    group_start_indices
                )[unordering_indices] / X_scale)
            }, simplify = FALSE)
        }
    }

    coefficient_log_likelihood <- vapply(coefficients_list, function(coeff) {
        return(ll_func(coeff[1:(p + 1)], coeff[(p + 2): (p + q + 2)], vect_y, mat_X, mat_Z, dist, a))
    }, numeric(1))

    # get number of model parameters (non-zero coefficients) and exclude intercepts - we have two intercepts
    # at the moment, this may be overstated => Zeng's approach of e^lambda may be better.
    df <- vapply(coefficients_list, function(x) { return(sum(x != 0) - 2) }, numeric(1))
    coefficients_bic <- vapply(seq_along(coefficient_log_likelihood), function(i) {
        return(df[i] * log(nrow(X)) - 2 * coefficient_log_likelihood[i])
    }, numeric(1))

    min_bic_index <- which.min(coefficients_bic)
    optimal_coefficients <- coefficients_list[[min_bic_index]]
    optimal_lambda <- lambda[min_bic_index]
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
                opt_lambda = optimal_lambda))
}

orthonormalise_full_matrix <- function(all_groups_matrix, group_start_indices) {
    num_cols <- ncol(all_groups_matrix)
    n <- nrow(all_groups_matrix)
    orthonormalised_matrix <- matrix(, nrow = n, ncol = num_cols)
    matrix_factors <- vector("list", length(group_start_indices))
    intercept_factors <- numeric(group_start_indices[2] - group_start_indices[1])
    for (i in group_start_indices[1]:(group_start_indices[2] - 1)) {
        orth_intercept <- svd_group_vector(all_groups_matrix[, i], n)
        orthonormalised_matrix[, i] <- orth_intercept[[1]]
        intercept_factors[i] <- orth_intercept[[2]]
    }
    matrix_factors[[1]] <- intercept_factors
    for (i in 2:(length(group_start_indices) - 1)) {
        # get the matrix:
        group_matrix <- all_groups_matrix[, group_start_indices[i]:(group_start_indices[i + 1] - 1), drop = FALSE]
        if (group_start_indices[i + 1] - 1 == group_start_indices[i]) {
            orth_group_matrix_and_factor <- svd_group_vector(group_matrix, n)
        } else {
            orth_group_matrix_and_factor <- svd_group_matrix(group_matrix, n)
        }
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

svd_group_vector <- function(vect, n) {
    svdresult <- svd(vect, nu = 0)
    vfactor <- sqrt(n) * svdresult$v * 1 / svdresult$d
    return(list(
        vect %*% vfactor,
        vfactor
    ))
}

l2_distance <- function(u, v) {
    return(sqrt(sum((u - v) ^ 2)))
}

group_coordinate_descent <- function(
    group_mat_orth, # orthonormal matrix of all groups/variables
    lambda,
    penalty_factors,
    initial_b,
    eps,
    group_start_indices,
    y,
    max_iterations = 1000,
    gamma = ifelse(penalty == "grSCAD", 4, 3),
    penalty = c("grLasso", "grMCP", "grSCAD", "grALasso")
) {
    # we need to output all of the coefficients from the different lambdae
    # we then calculate loss etc using bic later
    # 'coefficients' below may need to be a vector instead of a list
    coefficients_list <- vector("list", length(lambda))
    n <- nrow(group_mat_orth)
    penalty <- match.arg(penalty)

    if (penalty == "grLasso" || penalty == "grALasso") {
        threshold <- vector_soft_threshold
    } else if (penalty == "grMCP") {
        threshold <- vector_mcp_firm_threshold
    } else if (penalty == "grSCAD") {
        threshold <- vector_scad_firm_threshold
    }

    sdy <- sqrt(sum(y ^ 2) / nrow(group_mat_orth))
    tol <- eps * sdy
    # print(convergence_threshold)
    for (i in seq_along(lambda)) {
        # do not keep coefficients in group structure
        penalty_factors_i <- penalty_factors * lambda[i]
        r <- y - group_mat_orth %*% initial_b
        current_b <- initial_b
        for (iteration in 1:max_iterations) {
            # do the intercepts first. We treat them as separate, not as one group. We just 'store' them as one group.
            max_change <- 0
            lambda_0 <- penalty_factors_i[1]
            for (j in group_start_indices[1]:(group_start_indices[2] - 1)) {
                b <- current_b[j]
                X <- group_mat_orth[, j] # will not be changed
                z <- (1 / n) * t(X) %*% r + b
                new_b <- threshold(z, lambda_0, gamma)
                b_change <- new_b - b
                r <- r - X %*% (b_change) # is there an error in the paper by Breheny, which says to use t(X)?
                b <- new_b
                current_b[j] <- b
                max_change <- max(max_change, abs(b_change))
            }

            for (j in 2:(length(group_start_indices) - 1)) {
                lambda_j <- penalty_factors_i[j]
                b <- current_b[group_start_indices[j]:(group_start_indices[j + 1] - 1)]
                X <- group_mat_orth[, group_start_indices[j]:(group_start_indices[j + 1] - 1)] # will not be changed
                z <- (1 / n) * t(X) %*% r + b
                new_b <- threshold(z, lambda_j, gamma)
                b_change <- new_b - b
                r <- r - X %*% b_change # is there an error in the paper by Breheny, which says to use t(X)?
                b <- new_b
                current_b[group_start_indices[j]:(group_start_indices[j + 1] - 1)] <- b
                max_change <- max(max_change, max(abs(b_change)))
            }

            if (max_change < tol) {
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

create_lambda <- function(lambda_max, lambda_min, log.lambda, nlambda) {
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
    return(lambda)
}

mean_squared_error <- function(pseudo_X, coeff, actual_obs) {
    predictions <- pseudo_X %*% coeff
    return(sum((predictions - actual_obs)^2) / nrow(pseudo_X))
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

ll_func <- function(beta.count, beta.zero, y, X, Z, dist, a = 1)
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

    if (dist == "poisson") {
        ll <- try(sum(dZIP(y, mu = mu, sigma = pzero, log = TRUE)), silent = TRUE)
    } else if (dist == "negbin") {
        ll <- try(sum(dZINBI(y, mu = mu, sigma = 1 / a, nu = pzero, log = TRUE)), silent = TRUE)
    }

    # error-handling
    if (class(ll) == "try-error")
    {
        ll <- NA
    }
    return(ll)
}

# Taken from Lambert's paper
zip_log_likelihood <- function(X, Z, coeff, y) {
    x_coeff <- coeff[seq_len(ncol(X))]
    z_coeff <- coeff[(ncol(X) + 1):(ncol(X) + ncol(Z))]
    # indices of where y is zero
    indices_y0 <- which(y == 0)
    indices_yg0 <- which(y > 0)
    X_0 <- X[indices_y0, ]
    X_g0 <- X[indices_yg0, ]
    Z_0 <- Z[indices_y0, ]
    y_g0 <- y[indices_yg0]

    sum_log_y0 <- sum(log(exp(Z_0 %*% z_coeff) + exp(-exp(X_0 %*% x_coeff))))

    sum_all_y <- sum(log(1 + exp(Z %*% z_coeff)))

    sum_yg0 <- sum(y_g0 * (X_g0 %*% x_coeff) - exp(X_g0 %*% x_coeff))

    # sum over all y > 0 of (log(y!) = sum to y from 1 of logn)
    # we could possibly make this even faster
    log_sum <- cumsum(log(1:max(y_g0)))
    sum_yg0_factorial <- sum(log_sum[y_g0])

    return(sum_log_y0 + sum_yg0 - sum_all_y - sum_yg0_factorial)
}

zipll_deriv <- function(X, Z, coeff, y, k) {
    p1 <- ncol(X)
    x_coeff <- coeff[1:p1]
    z_coeff <- coeff[(p1 + 1):(p1 + ncol(Z))]
    # indices of where y is zero
    indices_y0 <- which(y == 0)
    indices_yg0 <- which(y > 0)
    X_0 <- X[indices_y0, ]
    X_g0 <- X[indices_yg0, ]
    Z_0 <- Z[indices_y0, ]
    y_g0 <- y[indices_yg0]

    # useful quantities
    A <- exp(X_0 %*% x_coeff)
    E <- exp(-A)
    J <- exp(Z_0 %*% z_coeff)

    A_g0 <- exp(X_g0 %*% x_coeff)
    J_all <- exp(Z %*% z_coeff)


    # count coefficient (beta) else zero coefficient (gamma)
    # we need to use k below
    if (k <= p1) {
        X_k_0 <- X[indices_y0, k]
        X_k_g0 <- X[indices_yg0, k]
        return(sum((-X_k_0 * A) / (1 + J / E)) + sum(X_k_g0 * (y_g0 - A_g0)))
    } else {
        Z_k_0 <- Z[indices_y0, k - p1]
        return(sum(Z_k_0 / (1 + E / J)) - sum(Z[, k - p1] / (1 + 1 / J_all)))
    }
}

zipll_dderiv <- function(X, Z, coeff, y, k) {
    p1 <- ncol(X)
    x_coeff <- coeff[1:p1]
    z_coeff <- coeff[(p1 + 1):(p1 + ncol(Z))]
    # indices of where y is zero
    indices_y0 <- which(y == 0)
    indices_yg0 <- which(y > 0)
    X_0 <- X[indices_y0, ]
    X_g0 <- X[indices_yg0, ]
    Z_0 <- Z[indices_y0, ]

    # useful quantities
    A <- exp(X_0 %*% x_coeff)
    E <- exp(-A)
    J <- exp(Z_0 %*% z_coeff)

    A_g0 <- exp(X_g0 %*% x_coeff)
    J_all <- exp(Z %*% z_coeff)

    # count coefficient (beta) else zero coefficient (gamma)
    # we need to use k below
    if (k <= p1) {
        X_k_0 <- X[indices_y0, k]
        X_k_g0 <- X[indices_yg0, k]
        return(sum((X_k_0^2 * E * A * (J * A + E + J)) / (E + J)^2) - sum(X_k_g0^2 * A_g0))
    } else {
        Z_k_0 <- Z[indices_y0, k - p1]
        return(sum(Z_k_0^2 * J * E / ((J + E)^2)) - sum(Z[, k - p1]^2 * J_all / ((1 + J_all)^2)))
    }
}

# when we do cycBAR there is no need to rearrange for groups

cycBAR <- function(X, Z, y, pseudo_X, pseudo_Y, X_scale, lambda, p, q, max_iterations, eps) {
    # no penalty factors here as this is individial variable selection so there is no need for group normalisation.
    n <- length(y)
    coefficients_list <- vector("list", length(lambda))
    sdy <- sqrt(sum(y ^ 2) / n)
    tol <- eps * sdy

    for (i in seq_along(lambda)) {
        lam <- lambda[i] / n
        # calculate ridge estimate
        coeff <- (solve(t(pseudo_X) %*% pseudo_X + lam * diag(ncol(pseudo_X))) %*% t(pseudo_X) %*% pseudo_Y) / X_scale
        for (iter in 1:max_iterations) {
            print(paste("iter", iter))
            new_coeff <- numeric(p + q + 2)
            # c2 <- hessian(function(coeff) {
            #     return(zip_log_likelihood(X, Z, coeff, y))
            # }, coeff)
            for (j in 1:(p + q + 2)) {
                c1j <- - zipll_deriv(X, Z, coeff, y, j)
                c2j <- - zipll_dderiv(X, Z, coeff, y, j)
                bj <- c2j * coeff[j] - c1j
                print(j)
                print(c2j)
                # print(c2[j, j])
                if (j != 1 && j != p + 2) {
                    if (abs(bj) < (2 * sqrt(c2j * lam))) {
                        new_coeff[j] <- 0
                    } else {
                        new_coeff[j] <- (bj + sign(bj) * sqrt(bj ^ 2 - 4 * c2j * lam)) / (2 * c2j)
                    }
                } else {
                    new_coeff[j] <- bj / c2j
                }
            }
            print(new_coeff)
            if (max(abs((new_coeff - coeff))) < tol) {
                break
            }
            coeff <- new_coeff
        }
        coefficients_list[[i]] <- coeff
    }

    return(coefficients_list)
}

maxgrad_bar <- function(X, Z, y, coefficients, p, q) {
    print("in maxgrad")
    print(ll_func(coefficients[1:(p + 1)], coefficients[(p + 2):(p + q + 2)], y, X, Z, dist = "poisson"))
    c1 <- - grad(function(coeff) {
        return(ll_func(coeff[1:(p + 1)], coeff[(p + 2):(p + q + 2)], y, X, Z, dist = "poisson"))
    }, coefficients)
    c2 <- - hessian(function(coeff) {
        return(ll_func(coeff[1:(p + 1)], coeff[(p + 2):(p + q + 2)], y, X, Z, dist = "poisson"))
    }, coefficients)
    print("c1")
    print(c1)
    print("c2")
    for (i in 1:(p + q + 2)) {
        print(c2[i, i])
    }
    lmax <- numeric(p + q + 2)
    for (i in 2:length(lmax)) {
        lmax[i] <- (c1[i]^2) / (4 * c2[i, i])
    }
    print("lmax:")
    print(lmax)
    print(which.max(lmax))
    return(max(lmax))
}
