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

// R, the third invariant of grad(u,v,w), which is just det(grad(u,v,w)).
// For 2D defined as: 
//      R = u_x * w_z - u_z * w_x
// For 3D defined as: 
//      R = u_x * (v_y * w_z - v_z * w_y) + u_y * (v_z * w_x - v_x * w_z) + u_z * (v_x * w_y - v_y * w_x)   
// This is used to find coherent structures correlated with eddies in turbulent flow
void R_invt(TArrayn::DTArray & R, TArrayn::DTArray & u,
        TArrayn::DTArray & v, TArrayn::DTArray & w, TArrayn::DTArray & temp1,
        TArrayn::DTArray & temp2, TArrayn::Grad * gradient_op,
        const string * grid_type, bool v_exist) {
    // Set-up
    S_EXP expan[3];
    assert(gradient_op);

    if (!v_exist) {
        //2D R
        // First term
        // Get u_x
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp1,false);
        // Get w_z
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp2,false);
        // u_x * w_z
        R = temp1 * temp2;
        
        // Second term
        // Get u_z
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp1,false);
        // Get w_x
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp2,false);
        // u_x * w_z
        R -= temp1 * temp2;
    }
    else {
        // 3D
        // First term
        // Get v_y
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp1,false);
        // Get w_z
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp2,false);
        // v_y * w_z
        temp1 = temp1 * temp2;
        // Get u_x
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp2,false);
        // v_y * w_z * u_x
        R = temp1 * temp2;
        // Get v_z
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp1,false);
        // Get w_y
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp2,false);
        // v_z * w_y
        temp1 = temp1 * temp2;
        // Get u_x
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp2,false);
        // - v_z * w_y * u_x
        R -= temp1 * temp2;

        // Second term
        // Get v_z
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp1,false);
        // Get w_x
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp2,false);
        // v_z * w_x      
        temp1 = temp1 * temp2;
        // Get u_y
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp2,false);
        // + v_z * w_x * u_y
        R += temp1 * temp2;
        // Get v_x
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp1,false);
        // Get w_z
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp2,false);
        // v_x * w_z
        temp1 = temp1 * temp2;
        // Get u_y
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp2,false);
        // - v_x * w_z * u_y
        R -= temp1 * temp2;

        // Third term
        // Get v_x
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp1,false);
        // Get w_y
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp2,false);
        // v_x * w_y
        temp1 = temp1 * temp2;
        // Get u_z
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp2,false);
        // + v_x * w_y * u_z
        R += temp1 * temp2;
        // Get v_y
        find_expansion(grid_type, expan, "v");
        gradient_op->setup_array(&v,expan[0],expan[1],expan[2]);
        gradient_op->get_dy(&temp1,false);
        // Get w_x
        find_expansion(grid_type, expan, "w");
        gradient_op->setup_array(&w,expan[0],expan[1],expan[2]);
        gradient_op->get_dx(&temp2,false);
        // v_y * w_x
        temp1 = temp1 * temp2;
        // Get u_z
        find_expansion(grid_type, expan, "u");
        gradient_op->setup_array(&u,expan[0],expan[1],expan[2]);
        gradient_op->get_dz(&temp2,false);
        // - v_y * w_x * u_z
        R -= temp1 * temp2;
    }
}
