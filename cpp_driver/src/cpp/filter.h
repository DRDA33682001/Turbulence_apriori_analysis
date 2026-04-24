/*
 * filter.h — Sharp spectral cutoff filter on a triply-periodic 3D grid.
 *
 * Physics reference
 *   Pope, Turbulent Flows (2000), Ch. 13:
 *     Eq. 13.22 defines the sharp spectral filter as the projection onto the
 *     subspace of Fourier modes with |k| <= kc = pi/delta.
 *     Eq. 13.93 (and the surrounding discussion) is the reminder that the
 *     SGS stress involves filter(ui*uj), NOT filter(ui)*filter(uj).
 *
 * Cross-check
 *   References/spectralDNS uses a wavenumber array K built from np.fft.fftfreq
 *   along the first two axes and np.fft.rfftfreq along the third (the r2c
 *   axis), then applies dealiasing / filtering by zeroing out coefficients
 *   whose |K| exceeds a threshold.  This implementation follows the same
 *   convention so a numpy-based a priori test would produce identical fields
 *   to round-off.
 *
 * Conventions
 *   - Domain is triply periodic with side length L = N * dx.
 *   - Physical wavenumbers are k_i = 2*pi * n_i / L, with n_i the signed
 *     integer wavenumber index.
 *   - FFTW r2c storage uses n_i in [0, N/2] on the contracted axis
 *     (Hermitian symmetry handles n_i < 0) and the standard signed layout
 *     n_i = i for i <= N/2, n_i = i - N for i > N/2 on the other axes.
 *   - Sharp cutoff keeps all modes with |k_phys| <= kc = pi/delta and zeros
 *     the rest.  The DC mode (n=0) is always retained.
 */

#pragma once

// Apply a sharp spectral cutoff filter at kc = pi/delta to a triply-periodic
// N x N x N field stored row-major in u_in.  Writes the filtered field into
// u_out.  u_in and u_out may alias safely (the transform uses internal
// scratch buffers).
//
//   u_in, u_out : flat 1D arrays of length N^3 (row-major, see config.h::idx)
//   N           : grid size per axis
//   delta       : filter width in the same physical units as dx
//   dx          : DNS grid spacing (domain length is L = N * dx)
//
// On the first call for a given N, this function builds FFTW plans with
// FFTW_MEASURE and loads/saves wisdom in results/fftw_wisdom.dat so that
// subsequent program invocations skip the measurement phase.
//
// Aborts the program if the filtered field has strictly more energy than
// the input (this would indicate a bug: sharp cutoff cannot add energy).
void spectral_filter_3d(const double* u_in,
                        double*       u_out,
                        int           N,
                        double        delta,
                        double        dx);

// Compute filter(ui * uj) correctly by first forming the pointwise product
// on the DNS grid and then applying the sharp spectral cutoff.  This is the
// "double-filtered" quantity that enters the SGS stress
//     tau_ij = filter(ui*uj) - filter(ui)*filter(uj)      [Pope Eq. 13.93]
// and it is NOT equal to filter(ui) * filter(uj).
//
//   ui, uj : input fields, length N^3 each
//   result : output, length N^3 (may alias ui or uj safely)
void filter_product(const double* ui,
                    const double* uj,
                    double*       result,
                    int           N,
                    double        delta,
                    double        dx);

// Apply a 3-D top-hat (physical-space box) filter on a triply-periodic grid.
// Output at each point is the unweighted mean of the input over a cube of
// edge length (2r+1) * dx centred on that point, with periodic wrap-around.
//
//   u_in, u_out : flat 1D arrays of length N^3 (row-major); may not alias
//   N           : grid size per axis
//   r           : box half-width in grid points (nominal box width = (2r+1)*dx)
//
// The top-hat filter is the classical LES test filter used in the Germano
// dynamic procedure (Germano et al. 1991, Lilly 1992).  It avoids the
// ill-posedness of the sharp spectral test filter at narrow alpha ratios,
// which produces wrong-signed mean Cs^2 in the Germano identity
// (Meneveau & Katz, ARFM 2000, sec. 3.2).
//
// Implementation: separable three-pass 1-D box filters using periodic
// prefix sums, O(N^3) total work independent of r.
void box_filter_3d(const double* u_in,
                   double*       u_out,
                   int           N,
                   int           r);

// Compute filter(ui * uj) where filter is the top-hat box filter.  Forms
// the pointwise product first, then applies the 3-D box filter.
//
//   ui, uj : input fields, length N^3 each
//   result : output, length N^3; may not alias ui or uj
void box_filter_product_3d(const double* ui,
                           const double* uj,
                           double*       result,
                           int           N,
                           int           r);
