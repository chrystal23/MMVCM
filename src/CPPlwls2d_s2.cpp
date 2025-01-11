#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect and sort

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd CPPlwls2d_s2( const double & bw, const std::string kernel_type, const Eigen::Map<Eigen::VectorXd> & win, 
                              const Eigen::Map<Eigen::MatrixXd> & tin, const Eigen::Map<Eigen::VectorXd> & yin, 
                              const Eigen::Map<Eigen::MatrixXd> & xin,
                              const Eigen::Map<Eigen::MatrixXd> & tout, const unsigned int & npoly = 1, 
                              const unsigned int & nder = 0){
  
  
  // Convenient constants
  const double invSqrt2pi=  1./(sqrt(2.*M_PI));
  const double factorials[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
  
  const unsigned int p = xin.cols();
  const unsigned int m = tin.cols();
  const unsigned int nTGrid = tin.rows();
  const unsigned int nUnknownPoints = tout.rows();
  Eigen::MatrixXd result(nUnknownPoints, p+1);
  
  // ========================
  // The checks start here:
  
  if(nTGrid == 0) {
    Rcpp::stop("The input T-grid has length zero.");
  }
  
  // Check that we have equal number of readings
  if( nTGrid != yin.size()){
    Rcpp::stop("The input Y-grid does not have the same number of points as input T-grid.");
  }
  
  if( nTGrid != win.size()){
    Rcpp::stop("The input weight vector does not have the same number of points as input T-grid.");
  }
  
  // Check that bandwidth is greater than zero
  if( bw <= 0.){
    Rcpp::stop("The bandwidth supplied for 1-D smoothing is not positive.");
  }
  
  // Check that the degree of polynomial used and the order of the derivative are reasonable
  if (npoly < nder){
    Rcpp::stop("The degree of polynomial supplied for 1-D smoothing is less than the order of derivative");
  }
  
  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels;
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4;
  possibleKernels["quar"]    = 5;
  
  // If the kernel_type key exists set KernelName appropriately
  int KernelName = 0;
  if ( possibleKernels.count( kernel_type ) != 0){
    KernelName = possibleKernels.find( kernel_type )->second; //Set kernel choice
  } else {
    // otherwise use "epan"as the kernel_type 
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;;
  }
  
  // Check that we do not have zero weights // Should do a try-catch here
  if ( !(win.all()) ){  // 
    Rcpp::warning("Cases with zero-valued windows maybe not be too safe.");
  }
  // Check if the first 5 elements are sorted // Very rough check // Will use issorted in the future
  if ( (tin(0,0)> tin(1,0)) ||  (tin(1,0)> tin(2,0)) ||  (tin(2,0)> tin(3,0)) ||  (tin(3,0)> tin(4,0)) ||  (tin(4,0)> tin(5,0)) ){
    Rcpp::stop("The T-grid used is not sorted. (or you have less than 6 points)");
  }
  
  // The checks end here.
  // ===================
  
  for (unsigned int i = 0; i != nUnknownPoints; ++i){
    //locating local window (LOL) (bad joke)
    std::vector <unsigned int> indx;
    const double* lower ;
    const double* upper ;
    
    //if the kernel is not Gaussian
    if ( KernelName != 3 && KernelName != 4) {
      //construct listX as vectors / size is unknown originally
      // for (unsigned int y = 0; y != nXGrid; ++y){ if ( std::fabs( xout(i) - xin(y) ) <= bw ) { indx.push_back(y); }  }
      // Get iterator pointing to the first element which is not less than tou(u)
      lower = std::lower_bound(tin.data(), tin.data() + nTGrid, tout(i, 0) - bw);
      upper = std::lower_bound(tin.data(), tin.data() + nTGrid, tout(i, 0) + bw);
      //  const unsigned int firstElement = lower - &tin[0];
      //  for (unsigned int xx1 =0; xx1 != upper-lower; ++xx1){
      //   indx.push_back(  firstElement+ xx1 );
      //  }
    } else {
      lower = tin.data();
      upper = tin.data() + nTGrid;
    }
    
    const unsigned int firstElement = lower - &tin(0, 0);
    for (unsigned int tt1 =0; tt1 != upper-lower; ++tt1){
      indx.push_back(  firstElement + tt1 );
    }
    
    // for(unsigned int r4=0; r4 != indx.size(); r4++){ Rcpp::Rcout << indx.at(r4) << ", "; } Rcpp::Rcout << "\n";
    
    unsigned int indxSize = indx.size();
    Eigen::VectorXd lw0(indxSize);
    Eigen::VectorXd ly0(indxSize);
    Eigen::MatrixXd lt0(indxSize, m);
    Eigen::MatrixXd lx0(indxSize, p);
    
    for (unsigned int y = 0; y != indxSize; ++y){
      lt0.row(y) = tin.row(indx[y]);
      lw0(y) = win(indx[y]);
      ly0(y) = yin(indx[y]);
      lx0.row(y) = xin.row(indx[y]);
    }

    std::vector <unsigned int> otherindx;
    if ( KernelName != 3 && KernelName != 4) {
      unsigned int flag = 1;
      for (unsigned int id = 0; id !=indxSize; ++id) {
        flag = 1;
        for (unsigned int l = 1; l != m; ++l) {
          if (lt0(id, l) < tout(i, l) - bw || lt0(id, l) > tout(i, l) + bw) {
            flag = 0;
            break;
          }
        }
        if (flag == 1) otherindx.push_back(id);
      }
    } else {
      for (unsigned int id = 0; id !=indxSize; ++id) otherindx.push_back(id);
    }

    indxSize = otherindx.size();
    Eigen::VectorXd temp(indxSize);
    Eigen::VectorXd lw(indxSize);
    Eigen::VectorXd ly(indxSize);
    Eigen::MatrixXd lt(indxSize, m);
    Eigen::MatrixXd lx(indxSize, p);
    
    for (unsigned int y = 0; y != indxSize; ++y){
      lt.row(y) = lt0.row(otherindx[y]);
      lw(y) = lw0(otherindx[y]);
      ly(y) = ly0(otherindx[y]);
      lx.row(y) = lx0.row(otherindx[y]);
    }
    
    
    Eigen::MatrixXd llt(indxSize, m);
    for (unsigned int l = 0; l != m; ++l){
        llt.col(l) = (lt.col(l).array() - tout(i, l)) * (1./bw);
    }
    
    //define the kernel used
    switch (KernelName){
      case 1: // Epan
        for (unsigned int y = 0; y != indxSize; ++y){
          temp(y) = 1;
          for (unsigned int l = 0; l != m; ++l) temp(y) *= (1- pow(llt(y, l), 2)) * 0.75 * (lw(y));
        }
        break;
      case 2 : // Rect
        for (unsigned int y = 0; y != indxSize; ++y){
          temp(y) = 1;
          for (unsigned int l = 0; l != m; ++l) temp(y) *= lw(y);
        }
        break;
      case 3 : // Gauss
        for (unsigned int y = 0; y != indxSize; ++y){
          temp(y) = 1;
          for (unsigned int l = 0; l != m; ++l) temp(y) *= exp(-.5*pow(llt(y, l), 2)) * invSqrt2pi *lw(y);
        }
        break;
      case 4 : // GausVar
        for (unsigned int y = 0; y != indxSize; ++y){
          temp(y) = 1;
          for (unsigned int l = 0; l != m; ++l) temp(y) *= exp(-.5*pow(llt(y, l), 2)) * invSqrt2pi *lw(y) * (1.25 - 0.25 * pow(llt(y, l), 2));
        }
        break;
      case 5 :  // Quar
        for (unsigned int y = 0; y != indxSize; ++y){
          temp(y) = 1;
          for (unsigned int l = 0; l != m; ++l) temp(y) *= pow((1.- pow(llt(y, l), 2)), 2) * (15./16.);
        }
        break;
    }

    
    if((p+1)*(m+1) >= indxSize){
      Rcpp::warning("Cannot do the estimation with less than (p+1)*(m+1)+1 points. tout id = %d", i);
      for (unsigned int y = 0; y != p+1; ++y) result(i, y) = std::numeric_limits<double>::quiet_NaN();
    } else {
      // make the design matrix
      Eigen::MatrixXd X(indxSize, (p+1)*(m+1));
      X.setOnes();
      for (unsigned int j = 1; j <= p; ++j){
        X.col(j) = lx.col(j-1).array();
      }
      for (unsigned int l = 0; l < m; ++l){
        X.col(p+1+l) = (tout(i, l) - lt.col(l).array()).array();
      }
      for (unsigned int j = 1; j <= p; ++j){
        for (unsigned int l = 0; l < m; ++l){
            X.col(p+1+j*m+l) = lx.col(j-1).array() * (tout(i, l) - lt.col(l).array()).array();
        }
      }
      
      Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() * X);
      // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
      Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
      
      //  Rcpp::Rcout << "lx: " << lx.transpose() << std::endl;
      //  Rcpp::Rcout << "ly: " << ly.transpose() << std::endl;
      //  Rcpp::Rcout << "temp: " << temp.transpose() << std::endl;
      //  Rcpp::Rcout << "llx: " << llx.transpose() << std::endl;
      //  Rcpp::Rcout << "xin: " << xin.transpose() << std::endl;
      //  Rcpp::Rcout << "yin: " << yin.transpose() << std::endl;
      //  Rcpp::Rcout << "xout: " << xout.transpose() << std::endl;
      //  Rcpp::Rcout << "X: " << X.transpose() << std::endl;
      //  Rcpp::Rcout << "beta: " << beta.transpose() << std::endl;
      //  Rcpp::Rcout << "factorials[nder]: " << factorials[nder]  << std::endl;
      //  Rcpp::Rcout << "pow (-1.0, nder): " << pow (-1.0, nder) << std::endl;
      
      for (unsigned int y = 0; y <= p; y++){
        result(i,y) = beta(y) * factorials[nder] *  std::pow(-1.0, int(nder));
      }
    }
  }
  return result;
}
