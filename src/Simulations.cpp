#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(".Simulation_CPP")]]
Eigen::VectorXd Simulation_CPP(Eigen::VectorXd eps, Eigen::MappedSparseMatrix<double> Wsparse, double alpha) {
  int n = eps.size();
  VectorXd eps_sq = eps.cwiseProduct(eps);

  SparseMatrix<double> I(n, n);
  I.setIdentity();

  SparseMatrix<double> M = I - Wsparse * eps_sq.asDiagonal();

  SparseLU<SparseMatrix<double>> solver;
  solver.compute(M);
  if(solver.info() != Success) {
    stop("Matrix M is singular.");
  }

  VectorXd alpha_vec = alpha * VectorXd::Ones(n);
  VectorXd h = solver.solve(alpha_vec);

  return h;
}

// [[Rcpp::export(".Simulation_CPP_LOG")]]
Eigen::VectorXd Simulation_CPP_LOG(Eigen::VectorXd eps, Eigen::MappedSparseMatrix<double> Wsparse, double alpha_0, double alpha_1) {
    int n = eps.size();

    VectorXd log_eps2 = eps.array().square().log().matrix();

    SparseMatrix<double> I(n, n);
    I.setIdentity();
    SparseMatrix<double> M = I - alpha_1 * Wsparse;

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(M);
    if(solver.info() != Success) {
        Rcpp::stop("Matrix (I - alpha_1 * W) is computationally singular in log-spARCH.");
    }

    VectorXd rhs = alpha_0 * VectorXd::Ones(n) + alpha_1 * Wsparse * log_eps2;
    VectorXd log_h = solver.solve(rhs);
    VectorXd Y = log_h.array().exp().sqrt().matrix().cwiseProduct(eps);
    return Y;
}

// [[Rcpp::export(".Simulation_CPP_E")]]
Eigen::VectorXd Simulation_CPP_E(Eigen::VectorXd eps, Eigen::MappedSparseMatrix<double> Wsparse, double alpha, double b) {
    int n = eps.size();
    VectorXd g_eps = (b * eps.array().abs().log()).matrix();
    VectorXd ln_h = alpha * VectorXd::Ones(n) + Wsparse * g_eps;
    VectorXd Y = ln_h.array().exp().sqrt().matrix().cwiseProduct(eps);
    return Y;
}



// [[Rcpp::export(".Simulation_spGARCH_CPP")]]
Eigen::VectorXd Simulation_spGARCH_CPP(double alpha, Eigen::VectorXd tau, Eigen::MappedSparseMatrix<double> W1, Eigen::MappedSparseMatrix<double> W2) {
  int n = tau.size();

  SparseMatrix<double> I(n, n);
  I.setIdentity();
  SparseMatrix<double> I_minus_W2 = I - W2;

  SparseLU<SparseMatrix<double>> solver;
  solver.compute(I_minus_W2);
  if(solver.info() != Success) {
    stop("Matrix I - W2 is computationally singular.");
  }

  VectorXd rhs = alpha * VectorXd::Ones(n) + W1 * tau;
  VectorXd X = solver.solve(rhs);

  return X;
}

// [[Rcpp::export(".d_h_d_eps_CPP")]]
NumericMatrix d_h_d_eps_CPP(const Eigen::VectorXd &eps,
                            const Eigen::VectorXd &h,
                            double alpha_scalar,
                            Eigen::MappedSparseMatrix<double> W1,
                            Eigen::MappedSparseMatrix<double> W2) {
  int n = eps.size();

  VectorXd eps2 = eps.cwiseProduct(eps);
  SparseMatrix<double> E(n,n);
  E.reserve(Eigen::VectorXi::Constant(n,1));
  for(int i=0;i<n;++i) E.insert(i,i) = eps2(i);

  SparseMatrix<double> I(n,n); I.setIdentity();
  SparseMatrix<double> M = I - W1 * E - W2;

  SparseLU<SparseMatrix<double>> lu;
  lu.compute(M);
  if(lu.info() != Success) Rcpp::stop("Matrix M is singular in d_h_d_eps_CPP");

  MatrixXd result = MatrixXd::Zero(n,n);
  VectorXd alpha_vec = VectorXd::Constant(n, alpha_scalar);
  VectorXd v = lu.solve(alpha_vec);

  SparseMatrix<double> S = I - W2;
  SparseLU<SparseMatrix<double>> luS;
  luS.compute(S);
  if(luS.info() != Success) Rcpp::stop("Matrix S is singular in d_h_d_eps_CPP");

  for(int j=0;j<n;++j){
    VectorXd e_j = VectorXd::Zero(n);
    e_j(j) = 1.0;
    VectorXd W1_col_j = W1 * e_j;
    VectorXd temp = lu.solve(W1_col_j);
    VectorXd W1_temp = W1 * temp;
    VectorXd w = luS.solve(W1_temp);
    double uj = v(j);
    result.col(j) = 2.0 * eps(j) * uj * w;
  }

  return wrap(result);
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
