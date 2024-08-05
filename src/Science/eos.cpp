#include "../Science.hpp"
#include "math.h"
#include "../Par_util.hpp"
#include "stdio.h"
#include "../Split_reader.hpp"
#include "../T_util.hpp"
#include "../Parformer.hpp"
#include "../Sorter.hpp"
#include <numeric>
#include <stdlib.h>

using blitz::Array;
using blitz::cos;
using namespace TArrayn;
using namespace NSIntegrator;
using namespace Transformer;


void eos(const string eos_type, TArrayn::DTArray & rho, TArrayn::DTArray & T, TArrayn::DTArray & S, double T0, double S0) {
    bool is_T0_valid = false;
    bool is_S0_valid = false;

    if (T0 >= -2) {
        is_T0_valid = true;
    }
    if (S0 >= 0) {
        is_S0_valid = true;
    }

    if (eos_type == "QUADEOS") {
        
        //For Debugging purposes
        /*
        if (master()) {
          	    fprintf(stdout,"You picked the QUADEOS option.\n");
                }
                MPI_Finalize(); exit(0);
        */
        //Call quadeos
        quadeos(rho,T);

    }
    else if (eos_type == "LINEOS") {
       if (is_T0_valid & is_S0_valid) {

        //For Debugging purposes
        /*
           if (master()) {
          	    fprintf(stdout,"You picked the LINEOS option.\n");
                }
                MPI_Finalize(); exit(0);
        */
           //Call lineos
        lineos(rho,T,S,T0,S0);
         }
       else {
           if (master()) {
          	fprintf(stderr,"Invalid option for background temperature or salinity. Exiting.\n");
           }
           MPI_Finalize(); exit(0);
       }
    }
    else if (eos_type == "NLEOS") {

        //For Debugging purposes
        /*
        if (master()) {
          	    fprintf(stdout,"You picked the NLEOS option.\n");
                }
                MPI_Finalize(); exit(0);
        */
        //Call nleos
        nleos(rho,T,S);

    }
    else {
       	    if (master()) {
			fprintf(stderr,"Invalid option received for eos type. Exiting.\n");
			}
                MPI_Finalize(); exit(0);
    }
}

 
