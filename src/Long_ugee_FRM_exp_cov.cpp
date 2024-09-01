#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP ugeeExp_interact_cov(Rcpp::List dij_all,
             Rcpp::List grp_all,
             Rcpp::List z_mat1,
             NumericVector wgroup,
             
             int t = 2,
             double tol = 1e-3,
             int maxiter = 50) {
  
  // This code is intended for t<=3.
  
  NumericVector g1 = grp_all[0];
  int Kx(unique(g1).size());
  int G_eff((Kx*(Kx+1)/2 - 1));
  
  int Kw(unique(wgroup).size());
  int W_eff((Kw*(Kw-1)/2));
  
  int p(1 + (t-1) + G_eff  + W_eff + 1);  // total number of parameters
  
  NumericMatrix v1 = dij_all[0];
  int n(v1.nrow());

  int iter(0);
  
  double err(1);
  
  VectorXd theta_00(VectorXd(p).setOnes());
  VectorXd theta(theta_00 * 0.7); // a smaller starting value
  VectorXd theta_new(theta * 1.0);
  
  MatrixXd delta(MatrixXd(t, p).setZero());
  VectorXd U(VectorXd(p).setZero());
  MatrixXd U_diff(MatrixXd(p, p).setZero());
  VectorXd fi_all(VectorXd(t).setZero());
  VectorXd h(VectorXd(t).setZero());
  
  VectorXd Xb(VectorXd(t).setZero());
  MatrixXd V_diag_mat(MatrixXd(t, t).setZero());

  while (err >= tol) {
    U.setZero();
    U_diff.setZero();
    
    for (int i = 0; i < n-1; ++i) {
      //for (int j = i+1; j < n; ++j) {
        for (int j = i+1; j < n; ++j) {
        
        fi_all.setZero();
        delta.setZero();
        
        
        for (int tt = 0; tt < t; ++tt) {
          
          //## fi ##
          NumericMatrix y_mat_tt = dij_all[tt];
          fi_all[tt] = y_mat_tt(i,j);
          
          //## Di ##
          NumericVector grp_tt = grp_all[tt];
          
          ///position: 0
          delta(tt, 0) = 1; // intercept
          
          ///position: 1  // This code is intended for t=2.
          if( (tt+1) == t){
            delta(tt, (t-1)) = 1;
          } else {
            delta(tt, (t-1)) = 0;
          }// time indicator
          
          ///position: 3
          int itermediate3ij((t-1) + Kx+(2*Kx-grp_tt[i])*(grp_tt[i]-1)/2+grp_tt[j]-grp_tt[i]-1);
          int itermediate3ji((t-1) + Kx+(2*Kx-grp_tt[j])*(grp_tt[j]-1)/2+grp_tt[i]-grp_tt[j]-1);
          
          int itermediate2((t-1) + grp_tt[i]-1);
          
          if (grp_tt[i] < grp_tt[j]) {
            delta(tt, itermediate3ij) = 1;
          } else if (grp_tt[i] > grp_tt[j]) {
            delta(tt, itermediate3ji) = 1;
          } else if (grp_tt[i] != 1) {
            ///position: 2
            delta(tt, itermediate2) = 1;
          }
          
          // ///position: 5
          // if (grp_tt[i] < grp_tt[j]) {
          //   delta(tt, (t-1) + G_eff +  Kx+(2*Kx-grp_tt[i])*(grp_tt[i]-1)/2+grp_tt[j]-grp_tt[i]-1) = delta(tt, (t-1));
          // } else if (grp_tt[i] > grp_tt[j]) {
          //   delta(tt, (t-1) + G_eff +  Kx+(2*Kx-grp_tt[j])*(grp_tt[j]-1)/2+grp_tt[i]-grp_tt[j]-1) = delta(tt, (t-1));
          // } else if (grp_tt[i] != 1) {
          //   ///position: 4
          //   delta(tt, (t-1) + G_eff +  grp_tt[i]-1) = delta(tt, (t-1));
          // }
          // 
          // Cate. cov vector wgroup1 [Time invariant]
          ///position: 6
          if (wgroup[i] != wgroup[j]) {
            delta(tt, p-2) = 1;
          } else {
            delta(tt, p-2) = 0;
          }// time indicator
          
          // Cont. cov matrix z_mat1
          NumericMatrix z_mat1_tt = z_mat1[tt];
          ///position: 7
          delta(tt, p-1) = z_mat1_tt(i,j);
          
          
        }
        
        Xb = (delta * theta);                // (t*p)*(p*1)=(t*1)
        h = exp(Xb.array());                 // (t*p)*(p*1)=(t*1)
        V_diag_mat = h.asDiagonal();     //(t*t)
        
        U += delta.transpose() * (fi_all - h); // (p*t)*(t*1)=(p*1)
        U_diff += delta.transpose() * V_diag_mat * delta;   // (p*t)*(t*t)*(t*p)=(p*p)
      }
    } // end of double for() loops
    
    MatrixXd U_diff_pinv = U_diff.completeOrthogonalDecomposition().pseudoInverse();
    theta_new = theta + U_diff_pinv * U;
    err = (theta_new - theta).lpNorm<Infinity>();
    theta = theta_new * 1.0;
    ++iter;
    
    if (iter > maxiter) {
      Rcerr << "Not converge!\n";
      break;
    }
  } // end of while()
  
  // Sandwich part:
  MatrixXd B(U_diff / (n * (n - 1) / 2));
  
  // Sigma_theta part:
  VectorXd v_i(VectorXd(p).setZero());
  MatrixXd Sigma_U(MatrixXd(p, p).setZero());
  
  
  for (int i = 0; i < n; ++i) {
    v_i.setZero();
    
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      
      fi_all.setZero();
      delta.setZero();
      
      for (int tt = 0; tt < t; ++tt) {
        
        //## fi ##
        NumericMatrix y_mat_tt = dij_all[tt];
        fi_all[tt] = y_mat_tt(i,j);
        
        //## Di ##
        NumericVector grp_tt = grp_all[tt];
        
        ///position: 0
        delta(tt, 0) = 1; // intercept
        
        ///position: 1  // This code is intended for t=2.
        if( (tt+1) == t){
          delta(tt, (t-1)) = 1;
        } else {
          delta(tt, (t-1)) = 0;
        }// time indicator
        
        ///position: 3
        int itermediate3ij((t-1) + Kx+(2*Kx-grp_tt[i])*(grp_tt[i]-1)/2+grp_tt[j]-grp_tt[i]-1);
        int itermediate3ji((t-1) + Kx+(2*Kx-grp_tt[j])*(grp_tt[j]-1)/2+grp_tt[i]-grp_tt[j]-1);
        
        int itermediate2((t-1) + grp_tt[i]-1);
        
        if (grp_tt[i] < grp_tt[j]) {
          delta(tt, itermediate3ij) = 1;
        } else if (grp_tt[i] > grp_tt[j]) {
          delta(tt, itermediate3ji) = 1;
        } else if (grp_tt[i] != 1) {
          ///position: 2
          delta(tt, itermediate2) = 1;
        }
        
        // ///position: 5
        // if (grp_tt[i] < grp_tt[j]) {
        //   delta(tt, (t-1) + G_eff +  Kx+(2*Kx-grp_tt[i])*(grp_tt[i]-1)/2+grp_tt[j]-grp_tt[i]-1) = delta(tt, (t-1));
        // } else if (grp_tt[i] > grp_tt[j]) {
        //   delta(tt, (t-1) + G_eff +  Kx+(2*Kx-grp_tt[j])*(grp_tt[j]-1)/2+grp_tt[i]-grp_tt[j]-1) = delta(tt, (t-1));
        // } else if (grp_tt[i] != 1) {
        //   ///position: 4
        //   delta(tt, (t-1) + G_eff +  grp_tt[i]-1) = delta(tt, (t-1));
        // }
        // 
        // Cate. cov vector wgroup1 [Time invariant]
        ///position: 6
        if (wgroup[i] != wgroup[j]) {
          delta(tt, p-2) = 1;
        } else {
          delta(tt, p-2) = 0;
        }// time indicator
        
        // Cont. cov matrix z_mat1
        NumericMatrix z_mat1_tt = z_mat1[tt];
        ///position: 7
        delta(tt, p-1) = z_mat1_tt(i,j);
        
      }
      
      Xb = (delta * theta);                // (t*p)*(p*1)=(t*1)
      h = exp(Xb.array());                 // (t*p)*(p*1)=(t*1)
      
      v_i += delta.transpose() * (fi_all - h) / (n - 1); // (p*t)*(t*1)=(p*1)
      
    }
    
    Sigma_U += 4 * v_i * v_i.adjoint() / (n - 1); // (p*1)*(1*p)=(p*p)
  } // end of double for() loops
  
  MatrixXd B_pinv = B.completeOrthogonalDecomposition().pseudoInverse();
  MatrixXd Sigma_theta(B_pinv * Sigma_U * B_pinv  / n);
  VectorXd sigma_theta(Sigma_theta.diagonal());
  
  return List::create(Named("theta") = wrap(theta),
                      Named("Sigma_theta") = wrap(Sigma_theta),
                      Named("sigma_theta") = wrap(sigma_theta));
}

