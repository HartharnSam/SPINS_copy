   /* Collection of various analysis routines that are general enough to be useful
   over more than one project */

#ifndef SCIENCE_HPP
#define SCIENCE_HPP 1

#include <blitz/array.h>
#include "TArray.hpp"
#include "NSIntegrator.hpp"


// Global arrays to store quadrature weights
extern Array<double,1> _quadw_x, _quadw_y, _quadw_z;

// Marek's Overturning Diagnostic
blitz::Array<double,3> overturning_2d(blitz::Array<double,3> const & rho,
      blitz::Array<double,1> const & zgrid, TArrayn::Dimension reduce = TArrayn::thirdDim );

// Read in a 2D file and interpret it as a 2D slice of a 3D array, for
// initialization with read-in-data from a program like MATLAB
void read_2d_slice(blitz::Array<double,3> & fillme, const char * filename,
                  int Nx, int Ny);

void read_2d_restart(blitz::Array<double,3>& fillme, const char* filename,
                  int Nx, int Ny);

// Vorticity
void compute_vort_x(TArrayn::DTArray & vortx, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type);
void compute_vort_y(TArrayn::DTArray & vorty, TArrayn::DTArray & u, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type);
void compute_vort_z(TArrayn::DTArray & vortz, TArrayn::DTArray & u, TArrayn::DTArray & v,
        TArrayn::Grad * gradient_op, const string * grid_type);
void compute_vorticity(TArrayn::DTArray & vortx, TArrayn::DTArray & vorty, TArrayn::DTArray & vortz,
        TArrayn::DTArray & u, TArrayn::DTArray & v, TArrayn::DTArray & w,
        TArrayn::Grad * gradient_op, const string * grid_type);

// Baroclinic Vorticity
void compute_baroclinic_vort_x(TArrayn::DTArray & barovortx, TArrayn::DTArray & T,
        TArrayn::Grad * gradient_op, const string * grid_type);
void compute_baroclinic_vort_y(TArrayn::DTArray & barovorty, TArrayn::DTArray & T,
        TArrayn::Grad * gradient_op, const string * grid_type);
void compute_baroclinic_vort(TArrayn::DTArray & barovort, TArrayn::DTArray & temp,
        TArrayn::DTArray & T, TArrayn::Grad * gradient_op, const string * grid_type, bool v_exist);

// Enstrophy density
void enstrophy_density(TArrayn::DTArray & enst, TArrayn::DTArray & u, TArrayn::DTArray & v,
        TArrayn::DTArray & w, TArrayn::Grad * gradient_op, const string * grid_type,
        const int Nx, const int Ny, const int Nz);

// Viscous dissipation
void dissipation(TArrayn::DTArray & diss, TArrayn::DTArray & u, TArrayn::DTArray & v,
        TArrayn::DTArray & w, TArrayn::Grad * gradient_op, const string * grid_type,
        const int Nx, const int Ny, const int Nz, const double visco);

// Background Potential Energy (BPE)
void compute_Background_PE(double & BPE_tot, TArrayn::DTArray & rho, TArrayn::DTArray & quad3,
        int Nx, int Ny, int Nz, double Lx, double Ly, double Lz, double g, double rho_0, int iter,
        bool dimensional_rho = false, bool mapped = false, Array<double,1> hill = Array<double,1>());

// Internal energy converted to BPE
void compute_BPE_from_internal(double & phi_i, TArrayn::DTArray & rho,
        double kappa_rho, double rho_0, double g, int Nz, bool dimensional_rho = false);

// Quadrature weights
void compute_quadweights(int szx, int szy, int szz,
      double Lx, double Ly, double Lz,
      NSIntegrator::DIMTYPE DIM_X, NSIntegrator::DIMTYPE DIM_Y,
      NSIntegrator::DIMTYPE DIM_Z);

const blitz::Array<double,1> * get_quad_x();
const blitz::Array<double,1> * get_quad_y();
const blitz::Array<double,1> * get_quad_z();

// find which expansion to use based on field and boundary conditions
void find_expansion(const string * grid_type, Transformer::S_EXP * expan, string deriv_filename);

// switch trig function
Transformer::S_EXP swap_trig( Transformer::S_EXP the_exp );

// Bottom slope
void bottom_slope(TArrayn::DTArray & Hprime, TArrayn::DTArray & zgrid,
        TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const int Nx, const int Ny, const int Nz);

// Top stresses
void top_stress_x(TArrayn::DTArray & stress_x, TArrayn::DTArray & u,
        TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const int Nz, const double visco);
void top_stress_y(TArrayn::DTArray & stress_y, TArrayn::DTArray & v,
        TArrayn::DTArray & temp, TArrayn::Grad * gradient_op,
        const string * grid_type, const int Nz, const double visco);

// Bottom stresses
void bottom_stress_x(TArrayn::DTArray & stress_x, TArrayn::DTArray & Hprime,
        TArrayn::DTArray & u, TArrayn::DTArray & w, TArrayn::DTArray & temp,
        TArrayn::Grad * gradient_op, const string * grid_type, const bool mapped,
        const double visco);
void bottom_stress_y(TArrayn::DTArray & stress_y, TArrayn::DTArray & Hprime,
        TArrayn::DTArray & v, TArrayn::DTArray & temp,
        TArrayn::Grad * gradient_op, const string * grid_type, const bool mapped,
        const double visco);

// Vortex stretching/tilting
void vortex_stretch_x(TArrayn::DTArray & vort_stretch, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op, const string * grid_type);
void vortex_stretch_y(TArrayn::DTArray & vort_stretch, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op, const string * grid_type);
void vortex_stretch_z(TArrayn::DTArray & vort_stretch, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op, const string * grid_type);

// Enstrophy production via vortex stretching/tilting
void enstrophy_stretch_production(TArrayn::DTArray & enst_prod, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::DTArray & temp3, TArrayn::Grad * gradient_op,
        const string * grid_type);

// Q/Second invariant of grad(u,v,w)
void Q_invt(TArrayn::DTArray & Q, TArrayn::DTArray & u,
         TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
         TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
         const string * grid_type);

// R/Third invariant of grad(u,v,w)
void R_invt(TArrayn::DTArray & R, TArrayn::DTArray & u,
         TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
         TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
         const string * grid_type, bool v_exist);


inline double nleos_inline(double T, double S){
   // Returns the density (kg/m^3)
   // MacDougall et. al. 2003  (JAOT 20)
   //
   // This EOS is a rational function taken at pressure = 0
   //
   // Constants for Numerator

    // First find any regions of negative salinity and set them to zero for this operation.

    if (S < 0) {
        S = 0;
    };

    const double N0 = 9.99843699e+02;
    const double N1 = 7.35212840e+00;
    const double N2 = -5.45928211e-02;
    const double N3 = 3.98476704e-04;
    const double N4 = 2.96938239e+00;
    const double N5 = -7.23268813e-03;
    const double N6 = 2.12382341e-03;
    const double N7 = 1.04004591e-02;
    const double N8 = 1.03970529e-07;
    const double N9 = 5.18761880e-06;
    const double N10 = -3.24041825e-08;
    const double N11 = -1.23869360e-11;

    const int p = 0; //Someday we could include pressure effects, and whoever wanted to do that would change p to an input variable from the model

    double numerator = N0 + T*(N1 + T*(N2 + N3*T)) + S*(N4 + N5*T + N6*S) + p*(N7 + N8*T*T + N9*S + p*(N10 + N11*T*T));

    // Constants for denominator
    const double D0 = 1.00000000e+00;
    const double D1 = 7.28606739e-03;
    const double D2 = -4.60835542e-05;
    const double D3 = 3.68390573e-07;
    const double D4 = 1.80809186e-10;
    const double D5 = 2.14691708e-03;
    const double D6 = -9.27062484e-06;
    const double D7 = -1.78343643e-10;
    const double D8 = 4.76534122e-06;
    const double D9 = 1.63410736e-09;
    const double D10 = 5.30848875e-06;
    const double D11 = -3.03175128e-16;
    const double D12 = -1.27934137e-17;

    double denominator = D0 + T*(D1 + T*(D2 + T*(D3 + D4*T))) + S*(D5 + T*(D6 + D7*T*T) + sqrt(S)*(D8 + D9*T*T)) + p*(D10 + p*T*(D11*T*T + D12*p));

   return numerator/denominator;
}
BZ_DECLARE_FUNCTION2(nleos_inline)

void nleos(TArrayn::DTArray & rho, TArrayn::DTArray & T,
         TArrayn::DTArray & S);

inline double compute_alpha(double T0, double S0){
    // Computes the thermal expansion coefficient at S0 and T0
    // Derivative is a finite difference for simplicity
    // Does not divide by rho_0, that is done in lineos()
    const double dT = 1e-08;
    double TpdT = T0 + dT;

    double alpha = (nleos_inline(TpdT,S0) - nleos_inline(T0,S0))/dT;

    return alpha;
}
BZ_DECLARE_FUNCTION(compute_alpha)

inline double compute_beta(double T0, double S0){
    // Computes the haline contraction coefficient at S0 and T0
    // Derivative is a finite difference for simplicity
    // Does not divide by rho_0, that is done in lineos()
    const double dS = 1e-08;
    double SpdS = S0 + dS;

    double beta = (nleos_inline(T0,SpdS) - nleos_inline(T0,S0))/dS;

    return beta;
}
BZ_DECLARE_FUNCTION(compute_beta)

inline double compute_rho0(double T0, double S0){
    // Computes the reference density at S0 and T0
    double rho_0 = nleos_inline(T0,S0);

    return rho_0;
}
BZ_DECLARE_FUNCTION(compute_rho0)
void lineos(TArrayn::DTArray & rho, TArrayn::DTArray & T,
         TArrayn::DTArray & S, const double & T0, const double & S0);

void quadeos(TArrayn::DTArray & rho, TArrayn::DTArray & T);


void eos(const string eos_type, TArrayn::DTArray & rho, TArrayn::DTArray & T, TArrayn::DTArray & S, double T0 = -10, double S0 = -2);

/*
inline double fresh_quad(double T){
   // Returns the density (kg/m^3) for water using simple quadratic fit to
   // MacDougall et. al. 2003  (JAOT 20)
   // Constants are determined via a quadratic fit preserving the Temperature of maximum density
   const double rho_max = 9.999744074665388e+02; // Density of freshwater at Tmd (deg C)
   //0 pressure relative to sea surface pressure.
   const double Tmd = 3.973973973973974; // Temperature of maximum density for freshwater at the surface
   const double C = -0.007641729398834; // Fitting Constant

   return rho_max + C*(T - Tmd)*(T - Tmd);
}
BZ_DECLARE_FUNCTION(fresh_quad)
*/

// lambda2, second eigenvalue of S^2+Omega^2
void compute_lambda2(TArrayn::DTArray & lambda2, TArrayn::DTArray & u,
    TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
    TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
    const string * grid_type, TArrayn::DTArray & A11, TArrayn::DTArray & A12,
    TArrayn::DTArray & A13, TArrayn::DTArray & A22, TArrayn::DTArray & A23,
    TArrayn::DTArray & A33);

/*  Creates and outputs data for plotting a bivariate histogram with NS*NT bins
 *  between tracers S1 and T1 to "filename.csv" (don't include file extension).
 */
struct QSPOptions {
  int NS;
  int NT;
  string filename;
  double S1_max;
  double S1_min;
  double T1_max;
  double T1_min;
  string T1_name;
  string S1_name;
};

struct QSPData {
  TArrayn::DTArray *u;
  TArrayn::DTArray *v;
  TArrayn::DTArray *w;
  TArrayn::DTArray *temp;
  TArrayn::DTArray *rho;
  TArrayn::DTArray *salinity;
  TArrayn::DTArray *custom_T1;
  TArrayn::DTArray *custom_S1;
  TArrayn::DTArray *xgrid;
  TArrayn::DTArray *ygrid;
  TArrayn::DTArray *zgrid;
  int Nx;
  int Ny;
  int Nz;
  int plotnum;
  bool mapped;
};

void QSPCount(QSPOptions qsp_options, QSPData qsp_data);

// Equation of state for seawater, polynomial fit from
// Brydon, Sun, Bleck (1999) (JGR)
//
inline double eqn_of_state(double T, double S){
   // Returns the density anomaly (kg/m^3) for water at the given
   // temperature T (degrees celsius) and salinity S (PSU?)

   // Constants are from table 4 of the above paper, at pressure 0
   // (This is appropriate since this is currently an incompressible
   // model)
   const double c1 = -9.20601e-2; // constant term
   const double c2 =  5.10768e-2; // T term
   const double c3 =  8.05999e-1; // S term
   const double c4 = -7.40849e-3; // T^2 term
   const double c5 = -3.01036e-3; // ST term
   const double c6 =  3.32267e-5; // T^3 term
   const double c7 =  3.21931e-5; // ST^2 term

   return c1 + c2*T + c3*S + c4*T*T + c5*S*T + c6*T*T*T + c7*S*T*T;
}

// Define a Blitz-friendly operator
BZ_DECLARE_FUNCTION2(eqn_of_state)

// Derivative of the equation of state WRT T

inline double eqn_of_state_dT(double T){
   // Returbs the derivative of the eqn_of_state function above WRT T
   // temperature T (degrees celsius), S assumed to be zero

   // Constants are from table 4 of the above paper, at pressure 0
   // (This is appropriate since this is currently an incompressible
   // model)

   // Kept the same names for these constants as in eqn_of_state
   const double c2 =  5.10768e-2; // T term, constant after derivative
   const double c4 = -7.40849e-3; // T^2 term, 2T after derivative
   const double c5 = -3.01036e-3; // ST term, S after derivative
   const double c6 =  3.32267e-5; // T^3 term, 3t^2 after derivative
   const double c7 =  3.21931e-5; // ST^2 term, 2ST after derivative

   return c2 + 2*c4*T + 3*c6*T*T;
}

// Define a Blitz-friendly operator
 BZ_DECLARE_FUNCTION(eqn_of_state_dT)

// Derivative of the equation of state WRT S

inline double eqn_of_state_dS(double T, double S){
   // Returns the density anomaly (kg/m^3) for water at the given
   // temperature T (degrees celsius) and salinity S (PSU?)

   // Constants are from table 4 of the above paper, at pressure 0
   // (This is appropriate since this is currently an incompressible
   // model)
   const double c3 =  8.05999e-1; // S term, constant after derivative
   const double c5 = -3.01036e-3; // ST term, T after derivative
   const double c7 =  3.21931e-5; // ST^2 term, T^2 after derivative

   return  c3 + c5*T + c7*T*T;
}

// Define a Blitz-friendly operator
BZ_DECLARE_FUNCTION2(eqn_of_state_dS)


inline double eqn_of_state_t(double T){
   // Specialize for freshwater with S=0
   return eqn_of_state(T,0.0);
}
BZ_DECLARE_FUNCTION(eqn_of_state_t)

#endif
