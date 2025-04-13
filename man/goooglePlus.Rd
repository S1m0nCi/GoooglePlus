\name{goooglePlus}
\alias{goooglePlus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
goooglePlus(data, xvars, zvars, yvar, group = 1:ncol(data), samegrp.overlap = TRUE, penalty = c("grLasso", "grMCP", "grSCAD", "grALasso"), dist = c("poisson", "negbin"), nlambda = 100, lambda, lambda_min = ifelse((nrow(data[, unique(c(xvars, zvars))]) > ncol(data[, unique(c(xvars, zvars))])), 1e-04, 0.05), lambda_max, log.lambda = TRUE, crit = "BIC", alpha = 1, eps = 0.001, max.iter = 1000, gmax = length(unique(group)), gamma = ifelse(penalty == "grSCAD", 4, 3), warn = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{xvars}{
%%     ~~Describe \code{xvars} here~~
}
  \item{zvars}{
%%     ~~Describe \code{zvars} here~~
}
  \item{yvar}{
%%     ~~Describe \code{yvar} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{samegrp.overlap}{
%%     ~~Describe \code{samegrp.overlap} here~~
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{lambda_min}{
%%     ~~Describe \code{lambda_min} here~~
}
  \item{lambda_max}{
%%     ~~Describe \code{lambda_max} here~~
}
  \item{log.lambda}{
%%     ~~Describe \code{log.lambda} here~~
}
  \item{crit}{
%%     ~~Describe \code{crit} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{gmax}{
%%     ~~Describe \code{gmax} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{warn}{
%%     ~~Describe \code{warn} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (data, xvars, zvars, yvar, group = 1:ncol(data), samegrp.overlap = TRUE, 
    penalty = c("grLasso", "grMCP", "grSCAD", "grALasso"), dist = c("poisson", 
        "negbin"), nlambda = 100, lambda, lambda_min = ifelse((nrow(data[, 
        unique(c(xvars, zvars))]) > ncol(data[, unique(c(xvars, 
        zvars))])), 1e-04, 0.05), lambda_max, log.lambda = TRUE, 
    crit = "BIC", alpha = 1, eps = 0.001, max.iter = 1000, gmax = length(unique(group)), 
    gamma = ifelse(penalty == "grSCAD", 4, 3), warn = TRUE) 
{
    y <- data.frame(data[, yvar])
    X <- data.frame(data[, xvars])
    pred.names <- union(xvars, zvars)
    names(X) <- paste(names(X), ".count", sep = "")
    names(y) <- yvar
    if (0 \%in\% group) {
        group_min <- min(group)
        group <- group + rep(group_min + 1, length(group))
    }
    group.x <- group[which(pred.names \%in\% xvars)]
    if (is.null(zvars)) {
        Z <- NULL
        group.z <- NULL
        data <- cbind.data.frame(y, X)
        xvars <- names(X)
        fit.formula <- as.formula(paste(yvar, "~", paste(paste(xvars, 
            collapse = "+"), "|1"), sep = ""))
    }
    else {
        Z <- data.frame(data[, zvars])
        names(Z) <- paste(names(Z), ".zero", sep = "")
        if (samegrp.overlap) {
            group.z <- group[which(pred.names \%in\% zvars)]
        }
        else {
            group.z <- max(group.x) + group[which(pred.names \%in\% 
                zvars)]
        }
        data <- cbind.data.frame(y, X, Z)
        xvars <- names(X)
        zvars <- names(Z)
        fit.formula <- as.formula(paste(yvar, "~", paste(paste(xvars, 
            collapse = "+"), "|", paste(zvars, collapse = "+")), 
            sep = ""))
    }
    p <- length(xvars)
    q <- length(zvars)
    fit.zero <- zeroinfl(fit.formula, dist = dist, data = data)
    a <- 1
    if (dist == "negbin") {
        a <- 1/fit.zero$theta
    }
    b2.mle <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)
    fit.coefficients <- c(fit.zero$coefficients$count, fit.zero$coefficients$zero)
    zeroinf_residuals <- fit.zero$residuals
    vcov <- fit.zero$vcov
    names(fit.coefficients)[1] <- "count.intercept"
    names(fit.coefficients)[p + 2] <- "zero.intercept"
    evcov <- eigen(vcov, symmetric = TRUE)
    inv_sqrt_diag <- diag(vapply(evcov$values, function(num) {
        num^(-1/2)
    }, numeric(1)), nrow = nrow(evcov$vectors), ncol = ncol(evcov$vectors))
    pseudo_factor <- evcov$vectors \%*\% inv_sqrt_diag \%*\% t(evcov$vectors)
    pseudo_X <- pseudo_factor
    pseudo_Y <- as.vector(pseudo_factor \%*\% b2.mle)
    X_std <- standardise(pseudo_X)
    pseudo_X <- X_std$scaled_X
    X_scale <- X_std$s
    fit.coefficients <- fit.coefficients * X_scale
    group <- c(0, group.x, 0, group.z)
    group_ordered_indices <- order(group)
    unordering_indices <- numeric(length(group))
    for (i in 1:length(group_ordered_indices)) {
        unordering_indices[group_ordered_indices[i]] <- i
    }
    pseudo_X <- pseudo_X[, group_ordered_indices]
    fit.coefficients <- fit.coefficients[group_ordered_indices]
    group <- sort(group)
    unique_groups <- unique(group)
    group_start_indices <- c(match(unique_groups, group), length(group) + 
        1)
    penalty_factors <- numeric(length(unique_groups))
    penalty_factors[1] <- 0
    for (i in 2:(length(group_start_indices) - 1)) {
        penalty_factors[i] <- sqrt(group_start_indices[i + 1] - 
            group_start_indices[i])
    }
    if (penalty == "grALASSO") {
        for (i in 2:(length(group_start_indices) - 1)) {
            group_norm <- sqrt(sum(fit.coefficients[group_start_indices[i]:group_start_indices[i + 
                1]]))
            if (group_norm == 0) {
                group_norm <- .Machine$double.eps
            }
            penalty_factors[i] <- penalty_factors[i]/group_norm
        }
    }
    group_matrices_orth_full <- orthonormalise_full_matrix(pseudo_X, 
        group_start_indices)
    group_matrix_orthonormal <- group_matrices_orth_full[[1]]
    orthonormalisation_factors <- group_matrices_orth_full[[2]]
    if (missing(lambda)) {
        if (missing(lambda_max)) {
            group_matrix_df <- as.data.frame(group_matrix_orthonormal)
            names(group_matrix_df) <- names(fit.coefficients)
            intercepts_fit <- glm(as.formula(paste("y ~ 0 + count.intercept + zero.intercept")), 
                data = cbind(data.frame(y = unlist(pseudo_Y)), 
                  group_matrix_df))
            zmax <- maxgrad(group_matrix_orthonormal, intercepts_fit$residuals, 
                group_start_indices, penalty_factors)/nrow(group_matrix_orthonormal)
            lambda_max <- zmax
        }
        if (log.lambda) {
            if (lambda_min == 0) {
                lambda <- c(exp(seq(log(lambda_max), log(0.001 * 
                  lambda_max), length = nlambda - 1)), 0)
            }
            else {
                lambda <- exp(seq(log(lambda_max), log(lambda_min * 
                  lambda_max), length = nlambda))
            }
        }
        else {
            if (lambda_min == 0) {
                lambda <- c(seq(lambda_max, 0.001 * lambda_max, 
                  length = nlambda - 1), 0)
            }
            else {
                lambda <- seq(lambda_max, lambda_min * lambda_max, 
                  length = nlambda)
            }
        }
    }
    orthonormalised_coefficients <- .Call("groupCoordinateDescent", 
        group_matrix_orthonormal, lambda, penalty_factors, fit.coefficients, 
        eps, group_start_indices, pseudo_Y, gamma, penalty, max.iter)
    coefficients_list <- lapply(orthonormalised_coefficients_list, 
        function(coeff) {
            return(deorthonormalise_coefficients(coeff, orthonormalisation_factors, 
                group_start_indices)[unordering_indices]/X_scale)
        })
    mat_X <- as.matrix(cbind(1, X))
    mat_Z <- as.matrix(cbind(1, Z))
    vect_y <- as.numeric(y[, 1])
    coefficient_log_likelihood <- vapply(coefficients_list, function(coeff) {
        return(ll.func(coeff[1:(p + 1)], coeff[(p + 2):(p + q + 
            2)], vect_y, mat_X, mat_Z, a))
    }, numeric(1))
    df <- vapply(coefficients_list, function(x) {
        return(sum(x != 0) - 2)
    }, numeric(1))
    coefficients_bic <- vapply(1:length(coefficient_log_likelihood), 
        function(i) {
            return(df[i] * log(nrow(X)) - 2 * coefficient_log_likelihood[i])
        }, numeric(1))
    min_bic_index <- which.min(coefficients_bic)
    optimal_coefficients <- coefficients_list[[min_bic_index]]
    optimal_lambda <- lambda[i]
    optimal_bic <- coefficients_bic[min_bic_index]
    count_coefficients <- optimal_coefficients[1:(p + 1)]
    zero_coefficients <- optimal_coefficients[(p + 2):(p + q + 
        2)]
    names(count_coefficients) <- names(fit.zero$coefficients$count)
    names(zero_coefficients) <- names(fit.zero$coefficients$zero)
    opt_params <- list(count = count_coefficients, zero = zero_coefficients)
    return(list(log_likelihood = coefficient_log_likelihood[min_bic_index], 
        bic = optimal_bic, coefficients = opt_params, opt_lambda = optimal_lambda))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
