/* Derivative case file for computing derivatives of existing fields */
// only able to return first and second derivatives

/* ------------------ Top matter --------------------- */

// Required headers
#include "../../BaseCase.hpp" // contains default class
#include "../../Options.hpp"  // config-file parser
#include "../../Science.hpp"  // Science content

#include <cassert>
#include <limits>
#include <string>

// Tensor variables for indexing
blitz::firstIndex ii;
blitz::secondIndex jj;
blitz::thirdIndex kk;

/* ------------------ Define parameters --------------------- */

// Grid scales
double Lx, Ly, Lz; // Grid lengths (m)
int Nx, Ny, Nz;    // Number of points in x, y, z
// Mapped grid?
bool mapped;
// Grid types
DIMTYPE intype_x, intype_y, intype_z;
string grid_type[3];
// Expansion types and indices
S_EXP expan[3];
const int x_ind = 0;
const int y_ind = 1;
const int z_ind = 2;

// QSP variables
double T1_max, S1_max, T1_min, S1_min;
std::string T1_name, S1_name;
int NT, NS;
std::string QSP_filename;
bool use_salinity, read_rho;
std::string salinity_filename, temp_filename, rho_filename;
std::string custom_T1_filename, custom_S1_filename;
const double double_max = std::numeric_limits<double>::max();

// current output number
int plotnum;

// Derivative options
int start_sequence; // output number to start taking derivatives at
int final_sequence; // output number to stop  taking derivatives at
int step_sequence;  // step between outputs to take derivatives
                    // streching/tilting?
bool v_exist;       // Does the v field exist?

/* ------------------ Adjust the class --------------------- */

class userControl : public BaseCase {
public:
  /* Initialize things */
  DTArray *xgrid, *ygrid, *zgrid; // Arrays for storing grid data

  /* Size of domain */
  double length_x() const { return Lx; }
  double length_y() const { return Ly; }
  double length_z() const { return Lz; }

  /* Resolution in X, Y, and Z */
  int size_x() const { return Nx; }
  int size_y() const { return Ny; }
  int size_z() const { return Nz; }

  /* Set expansions (FREE_SLIP, NO_SLIP (in vertical) or PERIODIC) */
  DIMTYPE type_x() const { return intype_x; }
  DIMTYPE type_y() const { return intype_y; }
  DIMTYPE type_z() const { return intype_z; }

  /* Set other things */
  int get_restart_sequence() const { return plotnum; }

  /* Read grid (if mapped) */
  bool is_mapped() const { return mapped; }
  void do_mapping(DTArray &xg, DTArray &yg, DTArray &zg) {
    init_grid_restart("x", "xgrid", xg);
    if (Ny > 1)
      init_grid_restart("y", "ygrid", yg);
    init_grid_restart("z", "zgrid", zg);
  }

  /* Read fields and do derivatives */
  void init_vels(DTArray &u, DTArray &v, DTArray &w) {
    // If user asks to grab the grids, we allocate arrays to store
    // them in memory
    if (mapped) {
      xgrid = alloc_array(Nx, Ny, Nz);
      if (Ny > 1) {
        ygrid = alloc_array(Nx, Ny, Nz);
      }
      zgrid = alloc_array(Nx, Ny, Nz);
      do_mapping(*xgrid, *ygrid, *zgrid);
    } else {
      // If user doesn't want to grab grids, we make sure not to
      // allocate arrays for them and to set the pointers to NULL.
      xgrid = NULL;
      ygrid = NULL;
      zgrid = NULL;
    }
    QSPOptions qsp_opts;
    qsp_opts.S1_name = S1_name;
    qsp_opts.T1_name = T1_name;
    qsp_opts.filename = QSP_filename;
    qsp_opts.NS = NS;
    qsp_opts.NT = NT;
    qsp_opts.S1_max = S1_max;
    qsp_opts.S1_min = S1_min;
    qsp_opts.T1_max = T1_max;
    qsp_opts.T1_min = T1_min;

    QSPData qsp_data;
    qsp_data.mapped = mapped;
    qsp_data.Nx = Nx;
    qsp_data.Ny = Ny;
    qsp_data.Nz = Nz;
    qsp_data.xgrid = xgrid;
    qsp_data.ygrid = ygrid;
    qsp_data.zgrid = zgrid;

    // Compute derivatives at each requested output
    for (plotnum = start_sequence; plotnum <= final_sequence;
         plotnum = plotnum + step_sequence) {

      // read in fields (if not already stored in memory)
      // Compute QSP data. The code promises to not mutate the arrays,
      // nor to make deep copies of them

      // u
      init_tracer_restart("u", u);
      // v
      if (v_exist) {
        init_tracer_restart("v", v);
      } else {
        if (master())
          fprintf(stdout, "No v field, setting v=0\n");
        v = 0;
      }
      // w
      init_tracer_restart("w", w);

      qsp_data.plotnum = plotnum;
      qsp_data.u = &u;
      qsp_data.v = &v;
      qsp_data.w = &w;
      qsp_data.temp = NULL;
      qsp_data.rho = NULL;
      qsp_data.salinity = NULL;
      qsp_data.custom_T1 = NULL;
      qsp_data.custom_S1 = NULL;

      DTArray *arr_salinity, *arr_temperature, *arr_rho, *custom_T1, *custom_S1;

      if (use_salinity) {
        arr_salinity = alloc_array(Nx, Ny, Nz);
        init_tracer_restart(salinity_filename, *arr_salinity);
        qsp_data.salinity = arr_salinity;
      }

      if (read_rho) {
        arr_rho = alloc_array(Nx, Ny, Nz);
        init_tracer_restart(rho_filename, *arr_rho);
        qsp_data.rho = arr_rho;
      }

      if (T1_name.compare("temp") == 0 || S1_name.compare("temp") == 0) {
        arr_temperature = alloc_array(Nx, Ny, Nz);
        init_tracer_restart(temp_filename, *arr_temperature);
        qsp_data.temp = arr_temperature;
      }

      if (T1_name.compare("custom") == 0) {
        custom_T1 = alloc_array(Nx, Ny, Nz);
        init_tracer_restart(custom_T1_filename, *custom_T1);
        qsp_data.custom_T1 = custom_T1;
      }

      if (S1_name.compare("custom") == 0) {
        custom_S1 = alloc_array(Nx, Ny, Nz);
        init_tracer_restart(custom_S1_filename, *custom_S1);
        qsp_data.custom_S1 = custom_S1;
      }

      QSPCount(qsp_opts, qsp_data);

      if (master()) {
        fprintf(stdout, "Completed QSP calculation for plotnum %d\n", plotnum);
      }
    }
  }

  // Constructor: Initialize local variables
  userControl() {
    compute_quadweights(size_x(), size_y(), size_z(), length_x(), length_y(),
                        length_z(), type_x(), type_y(), type_z());
  }
};

/* The ''main'' routine */
int main(int argc, char **argv) {
  /* Initialize MPI.  This is required even for single-processor runs,
     since the inner routines assume some degree of parallelization,
     even if it is trivial. */
  MPI_Init(&argc, &argv);

  /* ------------------ Define parameters from spins.conf ---------------------
   */
  options_init();

  option_category("Grid Options");
  add_option("Lx", &Lx, "Length of tank");
  add_option("Ly", &Ly, 1.0, "Width of tank");
  add_option("Lz", &Lz, "Height of tank");
  add_option("Nx", &Nx, "Number of points in X");
  add_option("Ny", &Ny, 1, "Number of points in Y");
  add_option("Nz", &Nz, "Number of points in Z");

  option_category("Grid mapping options");
  add_option("mapped_grid", &mapped, false, "Is the grid mapped?");

  string xgrid_type, ygrid_type, zgrid_type;
  add_option("type_x", &xgrid_type,
             "Grid type in X.  Valid values are:\n"
             "   FOURIER: Periodic\n"
             "   FREE_SLIP: Cosine expansion\n"
             "   NO_SLIP: Chebyhsev expansion");
  add_option("type_y", &ygrid_type, "FOURIER", "Grid type in Y");
  add_option("type_z", &zgrid_type, "Grid type in Z");

  option_category("QSP options");
  add_option("start_sequence", &start_sequence,
             "Sequence number to start taking derivatives at");
  add_option("final_sequence", &final_sequence,
             "Sequence number to stop  taking derivatives at");
  add_option("step_sequence", &step_sequence, 1,
             "Step between outputs to take derivatives");
  add_option("salinity_filename", &salinity_filename,
             "Base Filename of Salinity data (optional).");
  add_option("temp_filename", &temp_filename,
             "Base Filename of temperature data (optional).");
  add_option("T1_filename", &custom_T1_filename,
             "Base Filename of custom tracer for T1.");
  add_option("S1_filename", &custom_S1_filename,
             "Base Filename of custom tracer for S1.");
  add_option("T1", &T1_name, "u",
             "Name of tracer 1 for QSP. Valid values are rho,u,v,w,temp, "
             "salinity or ke");
  add_option("S1", &S1_name, "w",
             "Name of tracer 2 for QSP. Valid values are rho,u,v,w,temp, "
             "salinity or ke");
  add_option("T1_max", &T1_max, double_max,
             "Maximum explicit bin for T1 in QSP.");
  add_option("T1_min", &T1_min, -double_max,
             "Minimum explicit bin for T1 in QSP.");
  add_option("S1_max", &S1_max, double_max,
             "Maximum explicit bin for S1 in QSP.");
  add_option("S1_min", &S1_min, -double_max,
             "Minimum explicit bin for S1 in QSP.");
  add_option("salinity", &use_salinity, false,
             "Should salinity be read in from a file?.");
  add_option("read_rho", &read_rho, false,
             "Should rho be read in from a file?.");
  add_option("QSP_filename", &QSP_filename, "QSP_default",
             "Filename to save data to. Don't include file extension.");
  add_option("NS", &NS, 10, "Number of bins for tracer S");
  add_option("NT", &NT, 10, "Number of bins for tracer T");
  add_option("v_exist", &v_exist, "Does the v field exist?");
  // Parse the options from the command line and config file
  options_parse(argc, argv);

  /* ------------------ Adjust and check parameters --------------------- */
  /* Now, make sense of the options received.  Many of these values
     can be directly used, but the ones of string-type need further procesing.
   */

  // parse expansion types
  parse_boundary_conditions(xgrid_type, ygrid_type, zgrid_type, intype_x,
                            intype_y, intype_z);
  // vector of string types
  grid_type[x_ind] = xgrid_type;
  grid_type[y_ind] = ygrid_type;
  grid_type[z_ind] = zgrid_type;

  // adjust Ly for 2D
  if (Ny == 1 and Ly != 1.0) {
    Ly = 1.0;
    if (master())
      fprintf(stdout, "Simulation is 2 dimensional, "
                      "Ly has been changed to 1.0 for normalization.\n");
  }

  /* ------------------ Do stuff --------------------- */
  userControl mycode; // Create an instantiated object of the above class
  FluidEvolve<userControl> kevin_kh(&mycode);
  kevin_kh.initialize();
  MPI_Finalize();
  return 0;
}
