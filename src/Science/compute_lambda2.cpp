#include "../Science.hpp"
#include <blitz/array.h>
#include <string>
#include <math.h>

using namespace Transformer;
using namespace blitz;

// Calculates the second largest eigenvalue of S^2+Omega^2
// S and Omega are the symmetric and antisymmetric components of u_i,j, respectively
// Note: S^2+Omega^2 = 0.5*(u_{i,k}u_{k,j}+u_{k,i}u_{j,k})  
// This is used to find coherent structures correlated with eddies in turbulent flow
void compute_lambda2(TArrayn::DTArray & lambda2, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
        const string * grid_type, TArrayn::DTArray & A11, TArrayn::DTArray & A12,
        TArrayn::DTArray & A13, TArrayn::DTArray & A22, TArrayn::DTArray & A23,
        TArrayn::DTArray & A33) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);
    TArrayn::DTArray * A_ref; //This will hold references to the components
    TArrayn::DTArray * vel_mm; //This will hold the velocities
    TArrayn::DTArray * vel_nn; //This will hold the velocities
    std::string vel_labs[3] = {"u","v","w"};
    
    int ctr = 0;
    // Construct S^2+Omega^2
    for (int mm=0; mm<3; mm++){
        for (int nn=mm; nn<3; nn++){
            // Choose which component will be constructed
            if (ctr==0){A_ref = &A11;}
            else if (ctr==1){A_ref = &A12;}
            else if (ctr==2){A_ref = &A13;}
            else if (ctr==3){A_ref = &A22;}
            else if (ctr==4){A_ref = &A23;}
            else if (ctr==5){A_ref = &A33;}

            // A switch for the velocities
            // Useful for computing the u_{i,k}u_{j,k} terms
            
            if (mm==0){vel_mm = &u;}
            else if (mm==1){vel_mm = &v;}
            else if(mm==2){vel_mm = &w;}
            if (nn==0){vel_nn = &u;}
            else if (nn==1){vel_nn = &v;}
            else if (nn==2){vel_nn = &w;}
        
            // u_{i,k}u_{k,j}
            // Get u_{i,x}
            find_expansion(grid_type, expan, vel_labs[mm]);
            gradient_op->setup_array(vel_mm,expan[0],expan[1],expan[2]);
            gradient_op->get_dx(&temp1,false);
            // Get u_{x,j}
            find_expansion(grid_type, expan, "u");
            gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
            if (nn==0) { gradient_op->get_dx(&temp2,false); } 
            else if (nn==1) { gradient_op->get_dy(&temp2,false); }
            else if (nn==2) { gradient_op->get_dz(&temp2,false); }
            *A_ref = temp1 * temp2;
            // Get u_{i,y}
            find_expansion(grid_type, expan, vel_labs[mm]);
            gradient_op->setup_array(vel_mm,expan[0],expan[1],expan[2]);
            gradient_op->get_dy(&temp1,false);
            // Get u_{y,j}
            find_expansion(grid_type, expan, "v");
            gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
            if (nn==0) { gradient_op->get_dx(&temp2,false); } 
            else if (nn==1) { gradient_op->get_dy(&temp2,false); }
            else if (nn==2) { gradient_op->get_dz(&temp2,false); }
            *A_ref = temp1 * temp2;
            // Get u_{i,z}
            find_expansion(grid_type, expan, vel_labs[mm]);
            gradient_op->setup_array(vel_mm,expan[0],expan[1],expan[2]);
            gradient_op->get_dz(&temp1,false);
            // Get u_{z,j}
            find_expansion(grid_type, expan, "w");
            gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
            if (nn==0) { gradient_op->get_dx(&temp2,false); } 
            else if (nn==1) { gradient_op->get_dy(&temp2,false); }
            else if (nn==2) { gradient_op->get_dz(&temp2,false); }
            *A_ref = temp1 * temp2;
            
            
            // u_{k,i}u_{j,k}
            // Get u_{x,i}
            find_expansion(grid_type, expan, "u");
            gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
            if (mm==0) { gradient_op->get_dx(&temp1,false); } 
            else if (mm==1) { gradient_op->get_dy(&temp1,false); }
            else if (mm==2) { gradient_op->get_dz(&temp1,false); }
            // Get u_{j,x}
            find_expansion(grid_type, expan, vel_labs[nn]);
            gradient_op->setup_array(vel_nn,expan[0],expan[1],expan[2]);
            gradient_op->get_dx(&temp2,false);
            *A_ref = temp1 * temp2;
            // Get u_{y,i}
            find_expansion(grid_type, expan, "v");
            gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
            if (mm==0) { gradient_op->get_dx(&temp1,false); } 
            else if (mm==1) { gradient_op->get_dy(&temp1,false); }
            else if (mm==2) { gradient_op->get_dz(&temp1,false); }
            // Get u_{j,y}
            find_expansion(grid_type, expan, vel_labs[nn]);
            gradient_op->setup_array(vel_nn,expan[0],expan[1],expan[2]);
            gradient_op->get_dy(&temp2,false);
            *A_ref = temp1 * temp2;
            // Get u_{z,i}
            find_expansion(grid_type, expan, "w");
            gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
            if (mm==0) { gradient_op->get_dx(&temp1,false); } 
            else if (mm==1) { gradient_op->get_dy(&temp1,false); }
            else if (mm==2) { gradient_op->get_dz(&temp1,false); }
            // Get u_{j,z}
            find_expansion(grid_type, expan, vel_labs[nn]);
            gradient_op->setup_array(vel_nn,expan[0],expan[1],expan[2]);
            gradient_op->get_dz(&temp2,false);
            *A_ref = temp1 * temp2;
            
            *A_ref = (*A_ref)*0.5;
            ctr++;    
        }
    }
   
   // Now for the eigenvalue algorithm from: 
   // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
  
   // In the above article, we transform A into B=(1/p)*(A-qI) 
   // Then, lambda2 = q + 2*p*cos(arccos(det(B)/2)/3 + k*pi/3), k=0,1,2 
   // Define q
   temp1 = (A11+A22+A33)/3.0;
   // Define p 
   temp2 = pow2(A11-temp1) + pow2(A22-temp1) + pow2(A33-temp1) +2.0*(pow2(A12)+pow2(A13)+pow2(A23));
   temp2 = sqrt(temp2/6.0);

   // Transform A 
   A11 = (A11-temp1)/temp2;
   A22 = (A22-temp1)/temp2;
   A33 = (A33-temp1)/temp2;
   A12 = A12/temp2;
   A13 = A13/temp2;
   A23 = A23/temp2;
    
   // Calculate the determinant, using lambda2 as a dummy array
   lambda2 = 0.5*(A11*(A22*A33-A23*A23)-A12*(A12*A33-A13*A23)+A13*(A12*A23-A13*A22));

   // Since det(B)/2 feeds into acos, make sure it's between -1 and 1
   for (int i = lambda2.lbound(blitz::firstDim); i <= lambda2.ubound(blitz::firstDim); i++) { // outer loop over slowest-varying dimension
        for (int k = lambda2.lbound(blitz::thirdDim); k <= lambda2.ubound(blitz::thirdDim); k++) { // middle loop over next-slowest varying dimension
            for (int j = lambda2.lbound(blitz::secondDim); j <= lambda2.ubound(blitz::secondDim); j++) { // inner loop over fastest varying dimeion
                if (isnan(lambda2(i,j,k))) {lambda2(i,j,k) = 0;}
                else if (lambda2(i,j,k) < -1.0) {lambda2(i,j,k) = -1.0;}
                else if (lambda2(i,j,k) > 1.0) {lambda2(i,j,k) = 1.0;}
                // first conditional works as a safety check because the det is NaN iff temp2=0; in this case the 2nd eig is temp1
                  }}} 
   
   // Now calculate the actual lamba2. 
   lambda2 = temp1+2.0*temp2*cos(acos(lambda2)/3.0+4.0*M_PI/3.0);


}
