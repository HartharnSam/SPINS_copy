#include "../BaseCase.hpp"

// Bottom stresses
void BaseCase::point_ts(TArrayn::DTArray & rho, TArrayn::DTArray & xx,
        TArrayn::DTArray & zz, const std:vector<double> &z_locs, double x_locs, double time, int itercount, bool restarting) {

    // set-up
    blitz::firstIndex ii;
    blitz::secondIndex jj;
	blitz::thirdIndex kk;

	// Find the nearest index in x to the desired x position
    blitz::Array<double, 1> temp1 = blitz::abs(xx(ii) - x_loc);
    int x_loc_ii = blitz::minIndex(temp1);
	
	// Initialize a vector to store the nearest z indices
    std::vector<int> z_loc_ii(z_locs.size());

    // Find the nearest indices in z for the desired z positions
    for (size_t j = 0; j < z_locs.size(); j++) {
        blitz::Array<double, 1> temp2 = blitz::abs(zz(x_loc_ii, 0, kk) - z_locs[j]);
        z_loc_ii[j] = blitz::minIndex(temp2);
    }
	
	// Record the rho values at the specified x and z positions
    std::vector<double> pt(z_locs.size());
    for (size_t j = 0; j < z_locs.size(); j++) {
        pt[j] = rho(x_loc_ii, 0, z_loc_ii[j]);
    }

    // write to a stress diagnostic file
    if (master()) {
        FILE * point_ts_file = fopen("point_ts.txt","a");
        assert(point_ts_file);
        if ( itercount==0 and !restarting ) {
			fprintf(point_ts_file,"Time");
			for (size_t j = 0; j < z_locs.size(); j++) {
				fprintf(point_ts_file, ", pt%zu", j+1);
			}
			fprintf(point_ts_file, "\n");
		}
		fprintf(point_ts_file, "%.17f", time);
		for (size_t j = 0; j < z_locs.size(); j++) {
			fprintf(point_ts_file, ", pt%zu", j+1);
			}
		fprintf(point_ts_file, "\n");
		}
        fclose(point_ts_file);
    }
}

