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
   
// X-component of the baroclinic vorticity -g*d\rho/dy
void compute_baroclinic_vort_x(TArrayn::DTArray & barovortx, TArrayn::DTArray & T,
        TArrayn::Grad * gradient_op, const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // Setup for dT/dy
    find_expansion(grid_type, expan, "t");
    gradient_op->setup_array(&T,expan[0],expan[1],expan[2]);
    // get dT/dy
    gradient_op->get_dy(&barovortx,false);
    // use chainrule to calculate barovortx
    barovortx = (-1)*barovortx*eqn_of_state_dT(T);
}

