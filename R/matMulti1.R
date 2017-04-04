matMulti1 <-
function(P, G, C, X)  P %*% X %*% C %*% G %*% C %*% t(X) %*% P

#matMulti1 <- cxxfunction( signature(W_="numeric",X_="numeric",
#Y_="numeric", Z_="numeric"), 
#'arma::mat W = Rcpp::as<arma::mat>(W_);
#  arma::mat X = Rcpp::as<arma::mat>(X_);
#  arma::mat Y = Rcpp::as<arma::mat>(Y_);
#  arma::mat Z = Rcpp::as<arma::mat>(Z_);
 
#  return wrap(W*Z*Y * X * Y * trans(Z)*W);'
#  , plugin="RcppArmadillo")
