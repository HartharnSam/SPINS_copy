#include "../Science.hpp"
#include "boost/lexical_cast.hpp"
#include <algorithm>
#include <blitz/array.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <string>
#include <vector>

using namespace Transformer;
using namespace blitz;

class QSPVector {
private:
  double *data;
  int length;
  int rows;
  int cols;

public:
  explicit QSPVector(int length)
      : data((double *)calloc(length, sizeof(double))), length(length), rows(0),
        cols(0) {
    if (!data) {
      std::cout << "Error! Data could not be initialized.\n";
    }
  }
  QSPVector(int Nx, int Ny)
      : data((double *)calloc(Nx * Ny, sizeof(double))), length(Nx * Ny),
        rows(Nx), cols(Ny) {
    if (!data) {
      std::cout << "Error! Data could not be initialized.\n";
    }
  }
  double *raw() const { return data; }
  int Length() const { return length; }
  int Rows() const { return rows; }
  int Cols() const { return cols; }
  double operator[](int i) const {
    assert(i >= 0 && i < length);
    return data[i];
  }
  double operator()(int row, int col) const {
    assert(row >= 0 && row < rows);
    assert(col >= 0 && col < cols);
    return data[(row * cols) + col];
  }
  double &operator()(int row, int col) {
    assert(row >= 0 && row < rows);
    assert(col >= 0 && col < cols);
    return data[(row * cols) + col];
  }
  ~QSPVector() { free(data); }
};

enum QSPType {
  QSP_u,
  QSP_v,
  QSP_w,
  QSP_ke,
  QSP_temp,
  QSP_rho,
  QSP_salinity,
  QSP_custom,
};

void QSP_write(int local_rank, const QSPVector &local_hist,
               const QSPOptions &qsp_options, const QSPData &qsp_data) {
  if (local_rank == 0) {
    QSPVector glob_hist(qsp_options.NS, qsp_options.NT);
    MPI_Reduce(local_hist.raw(), glob_hist.raw(), // send and receive buffers
               qsp_options.NS * qsp_options.NT,   // Count
               MPI_DOUBLE,                        // datatype
               MPI_SUM, 0,      // Reduction operator and root process #
               MPI_COMM_WORLD); // Communicator
    string filename = qsp_options.filename + "." +
                      boost::lexical_cast<string>(qsp_data.plotnum) + ".csv";
    std::fstream outfile;
    outfile.open(filename.c_str(), std::ios_base::out);
    if (outfile.is_open()) {
      outfile << qsp_options.T1_max << ',' << qsp_options.T1_min << ','
              << qsp_options.S1_max << ',' << qsp_options.S1_min;
      for (int i = 4; i < qsp_options.NT; i++) {
        outfile << ',' << 0;
      }
      outfile << std::endl;
      for (int ii = 0; ii < qsp_options.NS; ii++) {
        outfile << glob_hist(ii, 0);
        for (int jj = 1; jj < qsp_options.NT; jj++) {
          outfile << ',' << glob_hist(ii, jj);
        }
        outfile << std::endl;
      }
    }
  } else {
    MPI_Reduce(local_hist.raw(), NULL,          // send and receive buffers
               qsp_options.NS * qsp_options.NT, // count
               MPI_DOUBLE,                      // datatype
               MPI_SUM, 0,      // Reduction operator and root process
               MPI_COMM_WORLD); // Communicator
  }
}

QSPType QSPConvert(const std::string &name) {
  QSPType converted_type;
  if (name.compare("u") == 0) {
    converted_type = QSP_u;
  } else if (name.compare("v") == 0) {
    converted_type = QSP_v;
  } else if (name.compare("w") == 0) {
    converted_type = QSP_w;
  } else if (name.compare("ke") == 0) {
    converted_type = QSP_ke;
  } else if (name.compare("temp") == 0) {
    converted_type = QSP_temp;
  } else if (name.compare("rho") == 0) {
    converted_type = QSP_rho;
  } else if (name.compare("salinity") == 0) {
    converted_type = QSP_salinity;
  } else if (name.compare("custom") == 0) {
    converted_type = QSP_custom;
  }
  return converted_type;
}

TArrayn::DTArray *QSPPtr(const QSPData &qsp_data, const QSPType &type) {
  TArrayn::DTArray *ptr = NULL;

  switch (type) {
  case QSP_u:
    ptr = qsp_data.u;
    break;
  case QSP_v:
    ptr = qsp_data.v;
    break;
  case QSP_w:
    ptr = qsp_data.w;
    break;
  case QSP_ke: // This is an odd case.
    break;
  case QSP_rho:
    ptr = qsp_data.rho;
    break;
  case QSP_temp:
    ptr = qsp_data.temp;
    break;
  case QSP_salinity:
    ptr = qsp_data.salinity;
    break;
  case QSP_custom: // Make sure to set the pointer yourself. Can't do it here.
    break;
  }

  return ptr;
}

// compute the minimum and maximum values for either tracer, then substitute
// any default values of +- infinity with the global min/max values instead.
void QSPMaxMin(const QSPType &T1_type, const QSPType &S1_type,
               TArrayn::DTArray *T1_ptr, TArrayn::DTArray *S1_ptr,
               QSPOptions &qsp_options, const QSPData &qsp_data, int i_low,
               int j_low, int k_low, int i_high, int j_high, int k_high) {
  double double_max = std::numeric_limits<double>::max();

  // If One of the variables is K.E or rho, we need to hand-roll the max / min
  // since, in general, the index of max(u) is the same as the index of max(v)
  if (T1_type == QSP_ke || T1_type == QSP_rho || S1_type == QSP_ke ||
      S1_type == QSP_rho) {
    double ke_max = -double_max, ke_min = double_max;
    double rho_max = -double_max, rho_min = double_max;
    // Main hand-rolled loop
    for (int i = i_low; i <= i_high; i++) {
      for (int j = j_low; j <= j_high; j++) {
        for (int k = k_low; k <= k_high; k++) {
          double tmp;
          if (T1_type == QSP_ke || S1_type == QSP_ke) {
            double ke_current = 0;
            tmp = (*qsp_data.u)(i, j, k);
            ke_current += tmp * tmp;
            if (qsp_data.Ny > 1) {
              tmp = (*qsp_data.v)(i, j, k);
              ke_current += tmp * tmp;
            }
            tmp = (*qsp_data.w)(i, j, k);
            ke_current += tmp * tmp;
            ke_current = 0.5 * ke_current;
            ke_max = ke_current > ke_max ? ke_current : ke_max;
            ke_min = ke_current < ke_min ? ke_current : ke_min;
          }
          if (T1_type == QSP_rho || S1_type == QSP_rho) {
            double rho_current;
            if (qsp_data.rho) {
              rho_current = (*qsp_data.rho)(i, j, k);
            } else if (qsp_data.temp && !qsp_data.salinity) {
              rho_current = eqn_of_state_t((*qsp_data.temp)(i, j, k));
            } else if (qsp_data.temp && qsp_data.salinity) {
              rho_current = eqn_of_state((*qsp_data.temp)(i, j, k),
                                         (*qsp_data.salinity)(i, j, k));
            } else {
              rho_current = 0;
            }
            rho_max = rho_current > rho_max ? rho_current : rho_max;
            rho_min = rho_current < rho_min ? rho_current : rho_min;
          }
        }
      }
    }
    double glob_ke_max, glob_ke_min, glob_rho_max, glob_rho_min;
    MPI_Allreduce(&ke_max, &glob_ke_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&ke_min, &glob_ke_min, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&rho_max, &glob_rho_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&rho_min, &glob_rho_min, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    switch (T1_type) {
    case QSP_ke:
      qsp_options.T1_max =
          qsp_options.T1_max == double_max ? glob_ke_max : qsp_options.T1_max;
      qsp_options.T1_min =
          qsp_options.T1_min == -double_max ? glob_ke_min : qsp_options.T1_min;
      break;
    case QSP_rho:
      qsp_options.T1_max =
          qsp_options.T1_max == double_max ? glob_rho_max : qsp_options.T1_max;
      qsp_options.T1_min =
          qsp_options.T1_min == -double_max ? glob_rho_min : qsp_options.T1_min;
      break;
    default:
      qsp_options.T1_max = qsp_options.T1_max == double_max
                               ? psmax(max(*T1_ptr))
                               : qsp_options.T1_max;
      qsp_options.T1_min = qsp_options.T1_min == -double_max
                               ? psmin(min(*T1_ptr))
                               : qsp_options.T1_min;
      break;
    }
    switch (S1_type) {
    case QSP_ke:
      qsp_options.S1_max =
          qsp_options.S1_max == double_max ? glob_ke_max : qsp_options.S1_max;
      qsp_options.S1_min =
          qsp_options.S1_min == -double_max ? glob_ke_min : qsp_options.S1_min;
      break;
    case QSP_rho:
      qsp_options.S1_max =
          qsp_options.S1_max == double_max ? glob_rho_max : qsp_options.S1_max;
      qsp_options.S1_min =
          qsp_options.S1_min == -double_max ? glob_rho_min : qsp_options.S1_min;
      break;
    default:
      qsp_options.S1_max = qsp_options.S1_max == double_max
                               ? psmax(max(*S1_ptr))
                               : qsp_options.S1_max;
      qsp_options.S1_min = qsp_options.S1_min == -double_max
                               ? psmin(min(*S1_ptr))
                               : qsp_options.S1_min;
      break;
    }
  } else { // !(cond1 || cond2) == !cond1 && !cond2
    qsp_options.S1_max = psmax(max(*S1_ptr));
    qsp_options.S1_min = psmin(min(*S1_ptr));
    qsp_options.T1_max = psmax(max(*T1_ptr));
    qsp_options.T1_min = psmin(min(*T1_ptr));
  }
}

void QSPCount(QSPOptions qsp_options, QSPData qsp_data) {

  int local_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

  // Find out what
  QSPType S1_type = QSPConvert(qsp_options.S1_name);
  QSPType T1_type = QSPConvert(qsp_options.T1_name);
  TArrayn::DTArray *S1_ptr = QSPPtr(qsp_data, S1_type);
  TArrayn::DTArray *T1_ptr = QSPPtr(qsp_data, T1_type);

  if (T1_type == QSP_custom) {
    T1_ptr = qsp_data.custom_T1;
  }

  if (S1_type == QSP_custom) {
    S1_ptr = qsp_data.custom_S1;
  }

  if ((!S1_ptr && S1_type != QSP_ke) || (!T1_ptr && T1_type != QSP_ke)) {
    std::cout << "Not enough data was provided for the requested tracer. "
                 "Aborting...\n";
    return;
  }

  int i_low, j_low, k_low, i_high, j_high, k_high;

  {
    TArrayn::DTArray *temp_ptr;
    if (S1_ptr) { // If S1 is not ke
      temp_ptr = S1_ptr;
    } else { // If S1 is ke we know u must exist
      temp_ptr = qsp_data.u;
    }
    i_low = temp_ptr->lbound(blitz::firstDim);
    j_low = temp_ptr->lbound(blitz::secondDim);
    k_low = temp_ptr->lbound(blitz::thirdDim);
    i_high = temp_ptr->ubound(blitz::firstDim);
    j_high = temp_ptr->ubound(blitz::secondDim);
    k_high = temp_ptr->ubound(blitz::thirdDim);
  }

  double double_max = std::numeric_limits<double>::max();
  if (qsp_options.T1_max == double_max || qsp_options.S1_max == double_max ||
      qsp_options.T1_min == -double_max || qsp_options.S1_min == -double_max) {
    QSPMaxMin(T1_type, S1_type, T1_ptr, S1_ptr, qsp_options, qsp_data, i_low,
              j_low, k_low, i_high, j_high, k_high);
  }

  double hS =
      (qsp_options.S1_max - qsp_options.S1_min) / (double)qsp_options.NS;
  double hT =
      (qsp_options.T1_max - qsp_options.T1_min) / (double)qsp_options.NT;
  double hS_inv = 1 / hS;
  double hT_inv = 1 / hT;

  QSPVector local_hist(qsp_options.NS, qsp_options.NT);
  QSPVector global_z_max(qsp_data.Nx, qsp_data.Ny);
  QSPVector global_z_min(qsp_data.Nx, qsp_data.Ny);
  // Find the range of Lz values per 2D-slice
  if (qsp_data.mapped) {
    QSPVector local_z_max(qsp_data.Nx, qsp_data.Ny);
    QSPVector local_z_min(qsp_data.Nx, qsp_data.Ny);
    //  We are slicing as if we are doing zgrid[i, j, :]
    for (int ii = i_low; ii <= i_high; ii++) {
      for (int jj = j_low; jj <= j_high; jj++) {
        // min is set to the highest possible value so it always gets changed
        double tmp_z_min = std::numeric_limits<double>::max();
        // max is set to the lowest possible value so it always gets changed
        double tmp_z_max = -tmp_z_min;
        double tmp;
        for (int kk = k_low; kk <= k_high; kk++) {
          tmp = (*qsp_data.zgrid)(ii, jj, kk);
          if (tmp > tmp_z_max) {
            tmp_z_max = tmp;
          } else if (tmp < tmp_z_min) {
            tmp_z_min = tmp;
          }
        }
        local_z_max(ii, jj) = tmp_z_max;
        local_z_min(ii, jj) = tmp_z_min;
      }
    }
    MPI_Allreduce(local_z_min.raw(), global_z_min.raw(),
                  qsp_data.Nx * qsp_data.Ny, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_z_max.raw(), global_z_max.raw(),
                  qsp_data.Nx * qsp_data.Ny, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  }

  // Main loop for QSP
  double Tval, Sval, tmp;
  for (int i = i_low; i <= i_high; i++) {
    for (int j = j_low; j <= j_high; j++) {
      for (int k = k_low; k <= k_high; k++) {

        switch (T1_type) {
        case QSP_ke:
          Tval = 0;
          tmp = (*qsp_data.u)(i, j, k);
          Tval += tmp * tmp;
          tmp = (*qsp_data.w)(i, j, k);
          Tval += tmp * tmp;
          if (qsp_data.Ny > 1) {
            tmp = (*qsp_data.v)(i, j, k);
            Tval += tmp * tmp;
          }
          Tval = 0.5 * Tval;
          break;
        case QSP_rho:
          tmp = (*T1_ptr)(i, j, k);
          Tval = eqn_of_state_t(tmp);
          break;
        default:
          Tval = (*T1_ptr)(i, j, k);
          break;
        }

        switch (S1_type) {
        case QSP_ke:
          Sval = 0;
          tmp = (*qsp_data.u)(i, j, k);
          Sval += tmp * tmp;
          tmp = (*qsp_data.w)(i, j, k);
          Sval += tmp * tmp;
          if (qsp_data.Ny > 1) {
            tmp = (*qsp_data.v)(i, j, k);
            Sval += tmp * tmp;
          }
          Sval = 0.5 * Sval;
          break;
        case QSP_rho:
          tmp = (*S1_ptr)(i, j, k);
          Sval = eqn_of_state_t(tmp);
          break;
        default:
          Sval = (*S1_ptr)(i, j, k);
          break;
        }

        int idxS = floor((Sval - qsp_options.S1_min) * hS_inv);
        int idxT = floor((Tval - qsp_options.T1_min) * hT_inv);
        idxS = std::max(std::min(idxS, qsp_options.NS - 1), 0);
        idxT = std::max(std::min(idxT, qsp_options.NT - 1), 0);

        double volume_weight;
        if (qsp_data.mapped) {
          // Calculate the Lz range
          double Lzmax_now = global_z_max(i, j);
          double Lzmin_now = global_z_min(i, j);

          // Calculate the arc length
          double arc, z_high, z_low, cos_high, cos_low;
          if (k > 0 && k < qsp_data.Nz - 1) {
            z_high = (double)(k + 1);
            z_low = (double)(k - 1);
          } else if (k == 0) {
            z_high = (double)(k + 1);
            z_low = (double)(k);
          } else if (k == qsp_data.Nz - 1) {
            z_high = (double)(k);
            z_low = (double)(k - 1);
          } else { // Failure
            std::cout << "k was out of bounds somehow. Failing...\n";
            return;
          }
          cos_high = std::cos(M_PI * z_high / (double)qsp_data.Nz);
          cos_low = std::cos(M_PI * z_low / (double)qsp_data.Nz);
          arc = 0.5 * (cos_low - cos_high);

          // Calculate the volume weight
          volume_weight = arc * (Lzmax_now - Lzmin_now);
        } else {
          volume_weight = 1.0;
        }

        local_hist(idxS, idxT) += volume_weight;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD); // Wait for everyone to finish
  QSP_write(local_rank, local_hist, qsp_options, qsp_data);
}
