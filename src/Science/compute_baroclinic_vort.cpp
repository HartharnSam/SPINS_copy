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
void compute_baroclinic_vort(TArrayn::DTArray & barovort, TArrayn::DTArray & temp,
        TArrayn::DTArray & T, TArrayn::Grad * gradient_op, const string * grid_type, bool v_exist ) {
    if ( v_exist ) {
       compute_baroclinic_vort_x(temp, T, gradient_op, grid_type);
       barovort = temp*temp; }
    else {
       barovort = 0;
    }
    compute_baroclinic_vort_y(temp, T, gradient_op, grid_type);
    barovort = barovort + temp*temp;
    barovort = sqrt(barovort);
}

