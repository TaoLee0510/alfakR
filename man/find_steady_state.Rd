% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction_functions.R
\name{find_steady_state}
\alias{find_steady_state}
\title{Find Steady State Karyotype Distribution}
\usage{
find_steady_state(lscape, p, Nmax = Inf)
}
\arguments{
\item{lscape}{data.frame. Must contain columns 'k' (character, karyotype string,
e.g., "2.2.1...") and 'mean' (numeric, fitness r_k). All karyotype strings
in \code{lscape$k} must be unique and represent karyotypes with the same number of chromosome types.}

\item{p}{numeric. Probability of single chromosome missegregation per division (0 <= p <= 1).}

\item{Nmax}{Optional limit to the number of missegregations allowable.}
}
\value{
A named numeric vector of steady-state frequencies, ordered according to
\code{lscape$k}. Returns \code{NULL} if calculation fails (e.g., invalid inputs,
eigen decomposition issues).
}
\description{
Calculates the steady-state frequency distribution of karyotypes by finding the
dominant eigenvector of the transition matrix M = t(W) \%*\% diag(r).
Assumes Rcpp functions for W matrix construction are compiled and available from the package.
}
\examples{
\dontrun{
# This example requires the C++ functions to be compiled and available
# within the alfakR package.

landscape_df <- data.frame(
  k = c("2.2", "3.1", "1.3"), # Example for 2 chromosome types
  mean = c(0.1, 0.12, 0.08)
)
prob_missegregation <- 0.001

# if (requireNamespace("alfakR", quietly = TRUE) &&
#     exists("rcpp_prepare_W_structure", where = "package:alfakR", mode="function")) {
#
#   ss_dist <- find_steady_state(lscape = landscape_df, p = prob_missegregation)
#   if (!is.null(ss_dist)) {
#     print("Steady-state distribution:")
#     print(ss_dist)
#     print(paste("Sum of frequencies:", sum(ss_dist)))
#   }
# } else {
#  message("Skipping find_steady_state example: alfakR not loaded or Rcpp fns not found.")
#  message("Run devtools::load_all() in your package project first, or install the package.")
# }
}
}
