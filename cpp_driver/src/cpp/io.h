/*
 * io.h — HDF5 I/O layer for the a priori LES SGS model evaluation.
 *
 * Provides:
 *   - VelocityField  : raw DNS velocity data (u1, u2, u3)
 *   - TensorField6   : symmetric 3x3 tensor stored as 6 components
 *   - ScalarField    : a single scalar quantity on the N^3 grid
 *
 * All arrays are flat 1D, row-major (C order), indexed by idx(i,j,k,N)
 * from config.h.  No physics lives here — pure I/O.
 *
 * HDF5 layout produced by this module:
 *
 *   Velocity input (read_velocity_hdf5)
 *     /               (root group)
 *       ATTR N        int64   — grid size
 *       ATTR dx       float64 — DNS grid spacing
 *       /u1           dataset (N,N,N) float64
 *       /u2           dataset (N,N,N) float64
 *       /u3           dataset (N,N,N) float64
 *
 *   Tensor output (write_tensor_hdf5)
 *     /<group_name>/
 *       ATTR N        int32
 *       /T_11  /T_12  /T_13  /T_22  /T_23  /T_33
 *           each dataset (N,N,N) float64
 *
 *   Scalar output (write_scalar_hdf5)
 *     /<dataset_name>   dataset (N,N,N) float64
 *
 *   Metrics output (write_metrics_hdf5)
 *     /<model_name>/
 *       ATTR correlation      float64
 *       ATTR backscatter_frac float64
 *       ATTR mean_nu_r        float64
 *       ATTR mean_Pi          float64
 *       ATTR mean_cos_theta   float64
 */

#pragma once

#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct VelocityField {
    std::vector<double> u1, u2, u3;  // flat 1D, row-major
    int    N;
    double dx;
};

struct TensorField6 {
    // Symmetric tensor stored as 6 components.
    // Order: 11, 12, 13, 22, 23, 33
    std::vector<std::vector<double>> components;  // 6 x N^3
    int         N;
    std::string name;
};

struct ScalarField {
    std::vector<double> data;  // flat 1D N^3
    int         N;
    std::string name;
};

// ---------------------------------------------------------------------------
// I/O functions
// ---------------------------------------------------------------------------

// Read velocity field from an HDF5 file produced by the JHTDB pipeline.
// Datasets 'u1', 'u2', 'u3' must be (N,N,N) float64 at the root group.
// Attributes 'N' (int64) and 'dx' (float64) are read from the root group.
// Prints shape and mean(u1) as a verification check.
VelocityField read_velocity_hdf5(const std::string& filepath);

// Write a symmetric tensor field (6 components, each N^3) into a named group.
// The file is created if it does not exist, otherwise it is opened for writing.
void write_tensor_hdf5(const TensorField6& T, const std::string& filepath,
                       const std::string& group_name);

// Write a scalar field (N^3) as a single (N,N,N) dataset.
// The file is created if it does not exist, otherwise it is opened for writing.
void write_scalar_hdf5(const ScalarField& s, const std::string& filepath,
                       const std::string& dataset_name);

// Write scalar SGS model diagnostics as attributes on a per-model group.
// Multiple models can be written to the same file under distinct group names.
void write_metrics_hdf5(const std::string& filepath,
                        const std::string& model_name,
                        double correlation,
                        double backscatter_frac,
                        double mean_nu_r,
                        double mean_Pi,
                        double mean_cos_theta);
