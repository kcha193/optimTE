matMulti <-
function(P, X, C) C %*% t(X) %*% P %*% X %*% C


#matMulti <- cxxfunction( signature(X="numeric",
#Y="numeric", Z="numeric"),  
#    ' arma::mat X = Rcpp::as<arma::mat>(X);
#      arma::mat Y = Rcpp::as<arma::mat>(Y);
#      arma::mat Z = Rcpp::as<arma::mat>(Z);
 
#      return wrap(Z*trans(Y) * X * Y * Z);', 
#      plugin="RcppArmadillo")
