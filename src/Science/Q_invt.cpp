#include "../Science.hpp"
/*#include "math.h"
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
*/
using namespace Transformer;

// Q, or the second invariant of grad(u,v,w)
// Defined as   u_x * v_y + v_y * w_z + u_x * w_z - u_y * v_x - v_z * w_y - u_z * w_x 
// This is used to find coherent structures correlated with eddies in turbulent flow
void Q_invt(TArrayn::DTArray & Q, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
        const string * grid_type) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    // First term
    // Get u_x
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&temp1,false);
    // Get v_y
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(&temp2,false);
    // u_x * v_y
    Q = temp1 * temp2;
    
    // Second term
    // Already have v_y, so get w_z
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(&temp1,false);
    // w_z * v_y
    Q += temp1 * temp2;

    // Third term
    // Already have w_z, so get u_z
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&temp2,false);
    // w_z * u_x
    Q += temp1 * temp2;
    
    // Fourth term
    // Get u_y
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(&temp1,false);
    // Get v_x
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&temp2,false);
    // - u_y * v_x
    Q -= temp1 * temp2;
    
    // Fifth term
    // Get v_z
    find_expansion(grid_type, expan, "v");
    gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(&temp1,false);
    // Get w_y
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dy(&temp2,false);
    // - u_y * v_x
    Q -= temp1 * temp2;
 
    // Sixth term
    // Get u_z
    find_expansion(grid_type, expan, "u");
    gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
    gradient_op->get_dz(&temp1,false);
    // Get v_x
    find_expansion(grid_type, expan, "w");
    gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
    gradient_op->get_dx(&temp2,false);
    // - u_z * w_x
    Q -= temp1 * temp2;
}
