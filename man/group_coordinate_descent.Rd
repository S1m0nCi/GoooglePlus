\name{group_coordinate_descent}
\alias{group_coordinate_descent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
group_coordinate_descent(group_mat_orth, lambda, penalty_factors, initial_b, initial_residuals, eps, group_start_indices, y, max_iterations = 1000, gamma = ifelse(penalty == "grSCAD", 4, 3), penalty = c("grLasso", "grMCP", "grSCAD", "grALasso"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group_mat_orth}{
%%     ~~Describe \code{group_mat_orth} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{penalty_factors}{
%%     ~~Describe \code{penalty_factors} here~~
}
  \item{initial_b}{
%%     ~~Describe \code{initial_b} here~~
}
  \item{initial_residuals}{
%%     ~~Describe \code{initial_residuals} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{group_start_indices}{
%%     ~~Describe \code{group_start_indices} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{max_iterations}{
%%     ~~Describe \code{max_iterations} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
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
function (group_mat_orth, lambda, penalty_factors, initial_b, 
    initial_residuals, eps, group_start_indices, y, max_iterations = 1000, 
    gamma = ifelse(penalty == "grSCAD", 4, 3), penalty = c("grLasso", 
        "grMCP", "grSCAD", "grALasso")) 
{
    coefficients_list <- vector("list", length(lambda))
    n <- nrow(group_mat_orth)
    penalty <- match.arg(penalty)
    if (penalty == "grLasso" || penalty == "grALasso") {
        threshold <- vector_soft_threshold
    }
    else if (penalty == "grMCP") {
        threshold <- vector_mcp_firm_threshold
    }
    else if (penalty == "grSCAD") {
        threshold <- vector_scad_firm_threshold
    }
    standard_deviation_y <- sqrt(sum(y^2)/nrow(group_mat_orth))
    convergence_threshold <- eps * standard_deviation_y
    for (i in 1:length(lambda)) {
        penalty_factors_i <- penalty_factors * lambda[i]
        r <- y - group_mat_orth \%*\% initial_b
        current_b <- initial_b
        for (iteration in 1:max_iterations) {
            max_change <- 0
            lambda_0 <- penalty_factors_i[1]
            for (j in group_start_indices[1]:(group_start_indices[2] - 
                1)) {
                b <- current_b[j]
                X <- group_mat_orth[, j]
                z <- (1/n) * t(X) \%*\% r + b
                new_b <- threshold(z, lambda_0, gamma)
                b_change <- new_b - b
                r <- r - X \%*\% (b_change)
                b <- new_b
                current_b[j] <- b
                max_change <- max(max_change, abs(b_change))
            }
            for (j in 2:(length(group_start_indices) - 1)) {
                lambda_j <- penalty_factors_i[j]
                b <- current_b[group_start_indices[j]:(group_start_indices[j + 
                  1] - 1)]
                X <- group_mat_orth[, group_start_indices[j]:(group_start_indices[j + 
                  1] - 1)]
                z <- (1/n) * t(X) \%*\% r + b
                new_b <- threshold(z, lambda_j, gamma)
                b_change <- new_b - b
                r <- r - X \%*\% b_change
                b <- new_b
                current_b[group_start_indices[j]:(group_start_indices[j + 
                  1] - 1)] <- b
                max_change <- max(max_change, max(abs(b_change)))
            }
            if (max_change < eps) {
                break
            }
        }
        coefficients_list[[i]] <- current_b
    }
    return(coefficients_list)
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
