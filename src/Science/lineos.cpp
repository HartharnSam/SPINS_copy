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

void lineos(TArrayn::DTArray & rho, TArrayn::DTArray & T, TArrayn::DTArray & S,
        const double & T0, const double & S0) { 

    double alpha = compute_alpha(T0,S0);
    double beta = compute_beta(T0,S0);
    double rho_0 = compute_rho0(T0,S0);

    rho = rho_0*(1 + alpha/rho_0*(T - T0) + beta/rho_0*(S - S0));

}
