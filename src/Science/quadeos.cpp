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


void quadeos(TArrayn::DTArray & rho, TArrayn::DTArray & T) {

   // Returns the density (kg/m^3) for water using simple quadratic fit to 
   // MacDougall et. al. 2003  (JAOT 20)
   // Constants are determined via a quadratic fit preserving the Temperature of maximum density
   const double rho_max = 9.999744074665388e+02; // Density of freshwater at Tmd (deg C)
   //0 pressure relative to sea surface pressure.
   const double Tmd = 3.973973973973974; // Temperature of maximum density for freshwater at the surface
   const double C = -0.007641729398834; // Fitting Constant

   rho = rho_max + C*(T - Tmd)*(T - Tmd);

}
