/*
 * sgs_exact.h — Exact (DNS-derived) sub-grid-scale stress and the strain-rate
 * tensor of the filtered velocity field.
 *
 * No modelling lives here. Everything in this file is derived directly from
 * the DNS velocity field by applying the sharp spectral filter and spectral
 * derivatives. The outputs are the "truth" against which SGS closures are
 * later compared in the a priori study.
 *
 * Physics references
 *   Pope, Turbulent Flows (2000):
 *     Eq. 13.93   — exact SGS stress  tau_ij = filter(u_i u_j) - u_bar_i u_bar_j
 *     p. 558      — the crucial note that filter(u_i u_j) != filter(u_i) * filter(u_j)
 *     Eq. 13.73   — rate-of-strain tensor S_ij = 0.5*(du_i/dx_j + du_j/dx_i)
 *     Eq. 13.74   — |S| = sqrt(2 S_ij S_ij)
 *
 * Cross-check
 *   references/spectralDNS computes spectral derivatives by multiplying the
 *   FFT coefficients by i*k_j (imaginary unit times the physical wavenumber)
 *   and inverse-transforming. This file follows exactly the same convention
 *   so a numpy-based verification would match to round-off.
 *
 * Notation
 *   vel, vel_bar : full DNS and filtered velocity fields (3 components each)
 *   tau          : exact SGS stress (symmetric, 6 unique components)
 *   tau_dev      : deviatoric part, tau_ij - (1/3) * trace(tau) * delta_ij
 *   g            : 9-component velocity gradient tensor, g[3*i+j] = du_bar_i/dx_j
 *   S            : symmetric strain-rate tensor of vel_bar (6 components)
 *   S_mag        : |S| = sqrt(2 S_ij S_ij), scalar field on the N^3 grid
 *
 * Symmetric tensor storage order (6 components): 11, 12, 13, 22, 23, 33
 * matching TensorField6 in io.h.
 */

#pragma once

#include "io.h"

#include <vector>

// ---------------------------------------------------------------------------
// Exact SGS stress tensor (Pope Eq. 13.93)
// ---------------------------------------------------------------------------
//
// Computes:
//   vel_bar_i = filter(u_i)                        for i = 1,2,3
//   tau_ij    = filter(u_i u_j) - vel_bar_i * vel_bar_j
//   tau_dev_ij = tau_ij - (1/3) * (tau_11 + tau_22 + tau_33) * delta_ij
//
// CRITICAL PHYSICS NOTE:
//   The order of operations matters. filter(u_i u_j) requires forming the
//   product on the FULL DNS grid and then filtering the result. This is
//   NOT the same as filter(u_i) * filter(u_j), which filters each factor
//   separately before multiplying. The difference between these two objects
//   is exactly the SGS stress tau_ij.
//   Reference: Pope (2000) p. 558, discussion immediately after Eq. 13.93.
//
// Inputs
//   vel   : DNS velocity field (u1, u2, u3) on the N^3 grid
//   delta : filter width (same physical units as dx)
//   dx    : DNS grid spacing
//
// Outputs (resized inside the function)
//   tau     : symmetric SGS stress, 6 unique components, each N^3
//   tau_dev : deviatoric SGS stress (same layout as tau)
//   vel_bar : filtered velocity field (3 components of length N^3)
//
void compute_exact_sgs_stress(const VelocityField& vel,
                              double               delta,
                              double               dx,
                              TensorField6&        tau,
                              TensorField6&        tau_dev,
                              VelocityField&       vel_bar);

// ---------------------------------------------------------------------------
// Velocity gradient tensor via spectral differentiation
// ---------------------------------------------------------------------------
//
// Computes the full 9-component velocity gradient tensor of the filtered
// field vel_bar:
//     g[3*i + j] = d(vel_bar_i) / d(x_j)      for i,j in {0,1,2}
//
// Convention: indices are 0-based, so i=0 -> u_bar_1, j=0 -> d/dx_1, etc.
//
// Spectral derivatives are EXACT on a periodic grid — the only error is
// round-off. They are obtained by forward FFT, multiplication by (i * k_j)
// with k_j the physical wavenumber, and inverse FFT. The Nyquist mode on
// each axis (if present in the r2c layout) is explicitly zeroed because
// its derivative would otherwise break Hermitian symmetry.
//
// Reference: references/spectralDNS uses the same i*k multiplication to
// compute derivatives on the triply-periodic box; see e.g. NS.py.
//
// Inputs
//   vel_bar : filtered velocity field
//   dx      : DNS grid spacing (domain length L = N * dx)
//
// Output (resized inside the function)
//   g : std::vector of 9 fields, each N^3, in row-major order as described
//
void compute_velocity_gradient_tensor(const VelocityField&              vel_bar,
                                      double                            dx,
                                      std::vector<std::vector<double>>& g);

// ---------------------------------------------------------------------------
// Strain-rate tensor from velocity gradient (Pope Eq. 13.73–13.74)
// ---------------------------------------------------------------------------
//
// S_ij = 0.5 * (g[3*i+j] + g[3*j+i])
// |S|  = sqrt(2 * S_ij * S_ij)    (full double sum over i, j)
//
// Symmetry of S is baked in by construction. The magnitude uses the double
// contraction 2 S_ij S_ij — the factor 2 matches Pope's definition of
// |S| that appears in the Smagorinsky eddy viscosity (Eq. 13.128).
//
// Inputs
//   g : 9-component velocity gradient tensor from
//       compute_velocity_gradient_tensor
//   N : grid size per axis (each g[k] has length N^3)
//
// Outputs (resized inside the function)
//   S     : symmetric strain rate, 6 unique components in TensorField6 order
//   S_mag : |S|, scalar field of length N^3 (non-negative by definition)
//
void compute_strain_rate_tensor(const std::vector<std::vector<double>>& g,
                                int                                     N,
                                TensorField6&                           S,
                                std::vector<double>&                    S_mag);
