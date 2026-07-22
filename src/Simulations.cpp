#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(".Simulation_CPP")]]
Eigen::VectorXd Simulation_CPP (Eigen::VectorXd eps,
                                  Eigen::MappedSparseMatrix<double> Wsparse,
                                  double alpha)
{
  int n = eps.size();

  VectorXd eps_sq = eps.array().square();
  VectorXd eta    = alpha * eps_sq;

  SparseMatrix<double> A = eps_sq.asDiagonal() * Wsparse;

  SparseMatrix<double> I_minus_A(n, n);
  I_minus_A.setIdentity();
  I_minus_A -= A;

  SparseLU<SparseMatrix<double>> solver;
  solver.compute(I_minus_A);
  VectorXd temp = solver.solve(eta);

  VectorXd h = alpha * VectorXd::Ones(n) + Wsparse * temp;

  return h;
}

// [[Rcpp::export(".Simulation_CPP_E")]]
Eigen::VectorXd Simulation_CPP_E (Eigen::ArrayXd eps, Eigen::MappedSparseMatrix<double> Wsparse, double alpha, double b)
{
  // ---------------- //
  //   Declarations   //
  // ---------------- //

  int n = eps.size();
  ArrayXd abs_eps_b(n);
  VectorXd g_eps(n);
  VectorXd Y(n);
  ArrayXd h(n);
  ArrayXd hexp(n);
  // ---------------- //
  //   Computations   //
  // ---------------- //

  abs_eps_b = pow(eps.abs(), b);
  g_eps = abs_eps_b.log();

  h = alpha * VectorXd::Ones(n) + Wsparse * g_eps;
  hexp = h.exp();

  Y = hexp.sqrt() * eps;

  return Y;
}



// [[Rcpp::export(".asdgCMatrix")]]
SEXP asdgCMatrix( SEXP IN ){
  // Søren Højsgaard
  // http://gallery.rcpp.org/articles/sparse-matrix-coercion/
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double
  MapMatd X(Rcpp::as<MapMatd>(IN));
  SpMat Xsparse = X.sparseView();              // Output: sparse matrix
  S4 Xout(wrap(Xsparse));                      // Output: as S4 object
  NumericMatrix Xin(IN);                       // Copy dimnames
  Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
  return(Xout);
}
