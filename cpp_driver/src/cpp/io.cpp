/*
 * io.cpp — HDF5 I/O implementation.
 *
 * Uses the HDF5 C API (libhdf5) directly — the CMake build links
 * HDF5_C_LIBRARIES, not the C++ wrapper.
 *
 * HDF5 C API reference:
 *   https://docs.hdfgroup.org/hdf5/v1_14/group___h5_f.html  (file)
 *   https://docs.hdfgroup.org/hdf5/v1_14/group___h5_d.html  (dataset)
 *   https://docs.hdfgroup.org/hdf5/v1_14/group___h5_a.html  (attribute)
 *   https://docs.hdfgroup.org/hdf5/v1_14/group___h5_g.html  (group)
 *   https://docs.hdfgroup.org/hdf5/v1_14/group___h5_s.html  (dataspace)
 *
 * Error handling policy:
 *   Every HDF5 call that returns an hid_t or herr_t is checked.
 *   On failure, all open HDF5 identifiers are closed before throwing.
 *   We suppress the HDF5 default error stack printer in open_or_create_file
 *   so that "file does not exist" probes do not pollute stderr.
 */

#include "io.h"
#include "config.h"   // idx(), N, DX, check_scalar, check_range

#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Silently test whether a file exists without triggering HDF5 error output.
static bool file_exists(const std::string& path) {
    struct stat st;
    return (stat(path.c_str(), &st) == 0);
}

// Open an HDF5 file for read/write; create it if it does not yet exist.
static hid_t open_or_create_file(const std::string& filepath) {
    hid_t fid;
    if (file_exists(filepath)) {
        fid = H5Fopen(filepath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (fid < 0)
            throw std::runtime_error("open_or_create_file: cannot open '" + filepath + "'");
    } else {
        fid = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fid < 0)
            throw std::runtime_error("open_or_create_file: cannot create '" + filepath + "'");
    }
    return fid;
}

// Write (or overwrite) a scalar float64 attribute on any open HDF5 object.
static void write_double_attr(hid_t loc_id, const char* name, double value) {
    // Delete pre-existing attribute so we can overwrite cleanly.
    if (H5Aexists(loc_id, name) > 0)
        H5Adelete(loc_id, name);

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr  = H5Acreate2(loc_id, name, H5T_IEEE_F64LE,
                              space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); throw std::runtime_error(std::string("write_double_attr: '") + name + "'"); }
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

// Write (or overwrite) a scalar int32 attribute.
static void write_int_attr(hid_t loc_id, const char* name, int value) {
    if (H5Aexists(loc_id, name) > 0)
        H5Adelete(loc_id, name);

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr  = H5Acreate2(loc_id, name, H5T_STD_I32LE,
                              space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); throw std::runtime_error(std::string("write_int_attr: '") + name + "'"); }
    H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

// Open or create a group; caller owns the returned hid_t.
static hid_t open_or_create_group(hid_t file_id, const std::string& name) {
    hid_t gid;
    if (H5Lexists(file_id, name.c_str(), H5P_DEFAULT) > 0) {
        gid = H5Gopen2(file_id, name.c_str(), H5P_DEFAULT);
    } else {
        gid = H5Gcreate2(file_id, name.c_str(),
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (gid < 0)
        throw std::runtime_error("open_or_create_group: cannot open/create '" + name + "'");
    return gid;
}

// Open or create a (N x N x N) float64 dataset; caller owns the returned hid_t.
static hid_t open_or_create_dataset_3d(hid_t loc_id, const char* dname, int N) {
    hid_t dset;
    if (H5Lexists(loc_id, dname, H5P_DEFAULT) > 0) {
        dset = H5Dopen2(loc_id, dname, H5P_DEFAULT);
    } else {
        hsize_t dims[3] = { (hsize_t)N, (hsize_t)N, (hsize_t)N };
        hid_t space = H5Screate_simple(3, dims, nullptr);
        dset = H5Dcreate2(loc_id, dname, H5T_IEEE_F64LE, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
    }
    if (dset < 0)
        throw std::runtime_error(std::string("open_or_create_dataset_3d: '") + dname + "'");
    return dset;
}

// Component names for the symmetric tensor in storage order: 11,12,13,22,23,33
static const char* const TENSOR_NAMES[6] = {
    "T_11", "T_12", "T_13", "T_22", "T_23", "T_33"
};

// ---------------------------------------------------------------------------
// read_velocity_hdf5
// ---------------------------------------------------------------------------

VelocityField read_velocity_hdf5(const std::string& filepath) {
    hid_t fid = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0)
        throw std::runtime_error("read_velocity_hdf5: cannot open '" + filepath + "'");

    VelocityField vf;

    // -- Attribute N (stored as int64 in file) --------------------------------
    {
        hid_t attr = H5Aopen(fid, "N", H5P_DEFAULT);
        if (attr < 0) { H5Fclose(fid); throw std::runtime_error("read_velocity_hdf5: attribute 'N' missing"); }
        long long N_ll = 0;
        H5Aread(attr, H5T_NATIVE_LLONG, &N_ll);
        H5Aclose(attr);
        vf.N = static_cast<int>(N_ll);
    }

    // -- Attribute dx (float64) -----------------------------------------------
    {
        hid_t attr = H5Aopen(fid, "dx", H5P_DEFAULT);
        if (attr < 0) { H5Fclose(fid); throw std::runtime_error("read_velocity_hdf5: attribute 'dx' missing"); }
        H5Aread(attr, H5T_NATIVE_DOUBLE, &vf.dx);
        H5Aclose(attr);
    }

    const long long total = (long long)vf.N * vf.N * vf.N;
    vf.u1.resize(total);
    vf.u2.resize(total);
    vf.u3.resize(total);

    // -- Datasets u1, u2, u3 --------------------------------------------------
    // The file stores each as (N,N,N) float64 in row-major order, which maps
    // directly onto our flat vectors without any stride adjustment.
    const char*          dnames[3] = { "u1",   "u2",   "u3"   };
    std::vector<double>* vecs[3]   = { &vf.u1, &vf.u2, &vf.u3 };

    for (int c = 0; c < 3; ++c) {
        hid_t dset = H5Dopen2(fid, dnames[c], H5P_DEFAULT);
        if (dset < 0) {
            H5Fclose(fid);
            throw std::runtime_error(std::string("read_velocity_hdf5: dataset '") + dnames[c] + "' not found");
        }
        herr_t err = H5Dread(dset, H5T_NATIVE_DOUBLE,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             vecs[c]->data());
        H5Dclose(dset);
        if (err < 0) {
            H5Fclose(fid);
            throw std::runtime_error(std::string("read_velocity_hdf5: H5Dread failed for '") + dnames[c] + "'");
        }
    }

    H5Fclose(fid);

    // -- Verification printout ------------------------------------------------
    // Mean of u1: quick sanity check. For isotropic turbulence the mean should
    // be near zero (mean flow removed in JHTDB dataset).
    double sum_u1 = 0.0;
    for (double v : vf.u1) sum_u1 += v;
    const double mean_u1 = sum_u1 / static_cast<double>(total);

    // RMS of each component
    double rms[3] = {0.0, 0.0, 0.0};
    std::vector<double>* rv[3] = { &vf.u1, &vf.u2, &vf.u3 };
    for (int c = 0; c < 3; ++c) {
        double sq = 0.0;
        for (double v : *rv[c]) sq += v * v;
        rms[c] = std::sqrt(sq / static_cast<double>(total));
    }

    std::printf("[io] read_velocity_hdf5: %s\n", filepath.c_str());
    std::printf("[io]   shape   : (%d, %d, %d)   N^3 = %lld\n",
                vf.N, vf.N, vf.N, total);
    std::printf("[io]   dx      : %.8f\n", vf.dx);
    std::printf("[io]   mean(u1): %+.6e\n", mean_u1);
    std::printf("[io]   rms(u1) : %.6e   rms(u2) : %.6e   rms(u3) : %.6e\n",
                rms[0], rms[1], rms[2]);

    // Mean-to-rms ratio: raw JHTDB sub-cube data retains a non-zero bulk mean,
    // so we allow up to 30 % before warning. Values above ~50 % suggest a
    // read error (wrong dataset, byte-swap, etc.).
    check_range("mean(u1) / rms(u1)", std::abs(mean_u1) / (rms[0] + 1e-14), 0.0, 0.30);

    return vf;
}

// ---------------------------------------------------------------------------
// write_tensor_hdf5
// ---------------------------------------------------------------------------

void write_tensor_hdf5(const TensorField6& T,
                       const std::string& filepath,
                       const std::string& group_name) {
    if (static_cast<int>(T.components.size()) != 6)
        throw std::runtime_error("write_tensor_hdf5: expected exactly 6 components");
    const long long expected = (long long)T.N * T.N * T.N;
    for (int c = 0; c < 6; ++c)
        if (static_cast<long long>(T.components[c].size()) != expected)
            throw std::runtime_error(std::string("write_tensor_hdf5: component ") +
                                     TENSOR_NAMES[c] + " has wrong length");

    hid_t fid = open_or_create_file(filepath);
    hid_t gid = open_or_create_group(fid, group_name);

    write_int_attr(gid, "N", T.N);

    for (int c = 0; c < 6; ++c) {
        hid_t dset = open_or_create_dataset_3d(gid, TENSOR_NAMES[c], T.N);
        herr_t err = H5Dwrite(dset, H5T_NATIVE_DOUBLE,
                              H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              T.components[c].data());
        H5Dclose(dset);
        if (err < 0) {
            H5Gclose(gid); H5Fclose(fid);
            throw std::runtime_error(std::string("write_tensor_hdf5: H5Dwrite failed for ") + TENSOR_NAMES[c]);
        }
    }

    H5Gclose(gid);
    H5Fclose(fid);

    std::printf("[io] write_tensor_hdf5: '%s' -> group '%s'  (%d^3 x 6 components)\n",
                filepath.c_str(), group_name.c_str(), T.N);
}

// ---------------------------------------------------------------------------
// write_scalar_hdf5
// ---------------------------------------------------------------------------

void write_scalar_hdf5(const ScalarField& s,
                       const std::string& filepath,
                       const std::string& dataset_name) {
    const long long expected = (long long)s.N * s.N * s.N;
    if (static_cast<long long>(s.data.size()) != expected)
        throw std::runtime_error("write_scalar_hdf5: data length mismatch");

    hid_t fid  = open_or_create_file(filepath);
    hid_t dset = open_or_create_dataset_3d(fid, dataset_name.c_str(), s.N);

    herr_t err = H5Dwrite(dset, H5T_NATIVE_DOUBLE,
                          H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          s.data.data());
    H5Dclose(dset);
    H5Fclose(fid);

    if (err < 0)
        throw std::runtime_error("write_scalar_hdf5: H5Dwrite failed for '" + dataset_name + "'");

    std::printf("[io] write_scalar_hdf5: '%s' -> dataset '%s'  (%d^3)\n",
                filepath.c_str(), dataset_name.c_str(), s.N);
}

// ---------------------------------------------------------------------------
// write_metrics_hdf5
// ---------------------------------------------------------------------------

void write_metrics_hdf5(const std::string& filepath,
                        const std::string& model_name,
                        double correlation,
                        double backscatter_frac,
                        double mean_nu_r,
                        double mean_Pi,
                        double mean_cos_theta) {
    hid_t fid = open_or_create_file(filepath);
    hid_t gid = open_or_create_group(fid, model_name);

    write_double_attr(gid, "correlation",      correlation);
    write_double_attr(gid, "backscatter_frac", backscatter_frac);
    write_double_attr(gid, "mean_nu_r",        mean_nu_r);
    write_double_attr(gid, "mean_Pi",          mean_Pi);
    write_double_attr(gid, "mean_cos_theta",   mean_cos_theta);

    H5Gclose(gid);
    H5Fclose(fid);

    std::printf("[io] write_metrics_hdf5: model '%s' -> '%s'\n",
                model_name.c_str(), filepath.c_str());
    std::printf("[io]   correlation=%.4f  backscatter_frac=%.4f  "
                "mean_nu_r=%+.4e  mean_Pi=%+.4e  mean_cos_theta=%.4f\n",
                correlation, backscatter_frac, mean_nu_r, mean_Pi, mean_cos_theta);
}
