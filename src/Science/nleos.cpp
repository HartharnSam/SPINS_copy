#include "../Science.hpp"
#include "math.h"
#include "../Par_util.hpp"
#include "stdio.h"
#include "../Split_reader.hpp"
#include "../T_util.hpp"
#include "../Parformer.hpp"
#include "../Sorter.hpp"
#include <numeric>

using blitz::Array;
using blitz::cos;
using namespace TArrayn;
using namespace NSIntegrator;
using namespace Transformer;


void nleos(TArrayn::DTArray & rho, TArrayn::DTArray & T, TArrayn::DTArray & S) {
   // Returns the density (kg/m^3) 
   // MacDougall et. al. 2003  (JAOT 20)
   //
   // This EOS is a rational function taken at pressure = 0
   //
   
    // First find any regions of negative salinity and set them to zero for this operation.
   

    //Numerator
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

    rho = N0 + T*(N1 + T*(N2 + N3*T)) + S*(N4 + N5*T + N6*S) + p*(N7 + N8*T*T + N9*S + p*(N10 + N11*T*T));

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
    
    rho = rho/(D0 + T*(D1 + T*(D2 + T*(D3 + D4*T))) + S*(D5 + T*(D6 + D7*T*T) + sqrt(max(S,0))*(D8 + D9*T*T)) + p*(D10 + p*T*(D11*T*T + D12*p)));
}

