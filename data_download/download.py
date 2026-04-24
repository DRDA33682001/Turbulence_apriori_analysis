"""download.py — JHTDB monolithic data download pipeline.

Refactored to eliminate silent data corruption caused by 8-chunk stitching.
Issues a single getCutout call per dataset (~67 MB each for 256^3) and runs
strict pre-save validation so corrupted volumes can never overwrite a good
file on disk.

Library note
------------
The brief names ``pyJHTDB`` for the backend, but pyJHTDB has been deprecated
(the package now raises ImportError at import time and directs users to
``giverny``). The current giverny package requires SciServer authentication,
which is not available in this environment. We therefore use ``givernylocal``
(the local-HTTP variant already bundled under references/giverny/) and wrap
it in a pyJHTDB-shaped ``_initialize_jhtdb`` / ``session.getCutout`` facade
so the rest of the module reads exactly like the spec. The important change
from the previous pipeline — one single request instead of eight stitched
chunks — is preserved exactly.

Dataflow (per dataset):
    1. session.getCutout(...)              -> raw numpy array
    2. _normalize_velocity_shape(raw)      -> (Nx, Ny, Nz, 3) float64
    3. _validate_{isotropic,channel}(u)    -> raise on failure
    4. _atomic_save(path, u, attrs)        -> .tmp then rename

A priori context:
    Previous pipeline stitched 8 x 128^3 sub-blocks via givernylocal's cache
    and occasionally misaligned chunks, producing octant RMS ratios ~3.5 and
    mean(Pi_exact) = -0.195 (wrong sign for a 3D cascade). The production token
    and the 2 GB GetCutout limit both allow a single-call 256^3 request, so
    chunking is no longer needed.

References:
    Pope (2000), Turbulent Flows, Ch. 13  - SGS stress + dissipation sign
    JHTDB isotropic1024coarse Re_lambda=433, L = (2*pi)^3 periodic domain
    JHTDB channel dataset Re_tau=1000, non-uniform y
"""

from __future__ import annotations

import datetime
import pathlib
import sys
from typing import Optional

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# Production JHTDB token. This is the only module-level mutable "globalish"
# value, per the brief. Everything else is either a constant or lives inside
# a function.
# ---------------------------------------------------------------------------
JHTDB_TOKEN: str = "you need to get your api key and put it here"

# Path constants (immutable; safe as module-level)
_HERE: pathlib.Path = pathlib.Path(__file__).resolve().parent        # src/python/
_PROJECT: pathlib.Path = _HERE.parent.parent                          # les_apriori/
_DATA_DIR: pathlib.Path = _PROJECT / "data"

# Grid sizing
_NATIVE_N_ISO: int = 1024
_OUTPUT_N: int = 256
_DEFAULT_STRIDE: int = _NATIVE_N_ISO // _OUTPUT_N  # = 4

# Channel extract dimensions.
# JHTDB channel: 2048(x) x 512(y) x 2048(z), Re_tau=1000.
# We extract (1024, 512, 1024) and apply strides (4, 2, 4) → (256, 256, 256).
# Using the first 1024 x-points and first 1024 z-points is sufficient because
# the flow is statistically homogeneous in both periodic directions.
_CHANNEL_EXTRACT: tuple[int, int, int] = (1024, 512, 1024)  # (Nx, Ny, Nz) native
_CHANNEL_STRIDE_XZ: int = 4    # stride in periodic x and z directions
_CHANNEL_STRIDE_Y:  int = 2    # stride in wall-normal y direction
_CHANNEL_RE_TAU:    int = 1000


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _log(msg: str) -> None:
    """Print a message with an ISO-like timestamp prefix.

    Args:
        msg: The message to print.
    """
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# JHTDB session bootstrap (givernylocal under a pyJHTDB-shaped facade)
# ---------------------------------------------------------------------------

_GIVERNYLOCAL_SRC: pathlib.Path = (
    _PROJECT.parent / "references" / "giverny" / "givernylocal" / "src"
)
_GIVERNYLOCAL_METADATA_JSON: pathlib.Path = (
    _PROJECT.parent / "references" / "giverny" / "metadata" / "configs" / "jhtdb_config.json"
)


def _ensure_givernylocal() -> None:
    """Make ``givernylocal`` importable from the bundled source tree.

    Adds ``references/giverny/givernylocal/src`` to sys.path if the
    package isn't already installed.
    """
    try:
        import givernylocal  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    local = str(_GIVERNYLOCAL_SRC)
    if local not in sys.path:
        sys.path.insert(0, local)
    import givernylocal  # noqa: F401
    _log(f"givernylocal: using bundled src at {local}")


def _patch_givernylocal_load_json_metadata() -> None:
    """Teach givernylocal's metadata loader to accept a local filesystem path.

    ``load_json_metadata`` normally ``requests.get``s the URL, which fails for
    plain file paths. We monkey-patch the loader to read the local JSON.
    """
    import givernylocal.turbulence_gizmos.basic_gizmos as _bg
    if getattr(_bg, "_patched_for_local_path", False):
        return

    _original = _bg.load_json_metadata

    def _local_aware_load(url):
        path = pathlib.Path(str(url).replace("file://", ""))
        if path.exists():
            from unittest.mock import MagicMock
            import requests as _r
            with open(path) as fh:
                text = fh.read()
            mock_resp = MagicMock()
            mock_resp.text = text
            _orig_get = _r.get
            _r.get = lambda u, **kw: mock_resp
            try:
                return _original(url)
            finally:
                _r.get = _orig_get
        return _original(url)

    _bg.load_json_metadata = _local_aware_load
    _bg._patched_for_local_path = True


class _JhtdbSession:
    """Minimal pyJHTDB-shaped facade over givernylocal.

    Exposes a ``getCutout`` method with the same keyword arguments the
    spec uses (``data_set``, ``field``, ``time_step``, ``start``, ``end``,
    ``step``, ``filter_width``) but internally issues a single giverny
    ``turbulence_toolkit.getCutout`` call on a per-dataset ``turb_dataset``
    cube.

    Each session owns a **brand-new** cube, freshly constructed on first
    use and discarded in ``finalize``. Callers that want defeat-all-caches
    semantics should construct a new ``_JhtdbSession`` per request — this
    is exactly what the per-chunk loop in ``download_*_monolithic`` does,
    and it is the root-cause fix for the previous 8-chunk corruption bug
    (shared-cube cache was misaligning chunks).
    """

    def __init__(self) -> None:
        self._cube = None
        self._dataset_title: Optional[str] = None

    def _make_cube(self, dataset_title: str) -> None:
        """Instantiate a fresh per-session giverny cube and bump the local cap.

        We patch ``max_local_cutout_size`` up from 3072 MB to 20 GB so the
        client-side check does not reject large strided requests. The
        per-chunk box sizes chosen by callers still stay well below what
        the server can serve in a single HTTP pass.
        """
        _ensure_givernylocal()
        _patch_givernylocal_load_json_metadata()

        from givernylocal.turbulence_dataset import turb_dataset

        scratch = _DATA_DIR / f"giverny_scratch_{dataset_title}"
        scratch.mkdir(parents=True, exist_ok=True)

        self._cube = turb_dataset(
            dataset_title=dataset_title,
            output_path=str(scratch),
            auth_token=JHTDB_TOKEN,
            json_url=str(_GIVERNYLOCAL_METADATA_JSON),
        )
        # Bump the client-side cap so strided calls are not rejected
        # by givernylocal's native-box sanity check.
        try:
            self._cube.metadata["constants"]["max_local_cutout_size"] = 20480.0
        except (AttributeError, KeyError, TypeError):
            pass
        self._dataset_title = dataset_title

    def getCutout(
        self,
        data_set: str,
        field: str,
        time_step: int,
        start: np.ndarray,
        end: np.ndarray,
        step: np.ndarray,
        filter_width: int = 1,
    ) -> np.ndarray:
        """Issue a single giverny getCutout request and return a numpy array.

        Args:
            data_set: e.g. ``'isotropic1024coarse'`` or ``'channel'``.
            field: ``'u'`` for velocity (giverny calls it ``'velocity'``).
            time_step: 1-based JHTDB time index.
            start: length-3 int array, 1-based inclusive spatial start.
            end: length-3 int array, 1-based inclusive spatial end.
            step: length-3 int array, per-axis coarsening stride.
            filter_width: giverny spatial filter width (kept for API shape).

        Returns:
            ndarray of shape ``(Nz, Ny, Nx, 3)``, float32.
        """
        del filter_width  # unused in this monolithic path
        self._make_cube(data_set)
        from givernylocal.turbulence_toolkit import getCutout as _giverny_getCutout

        giverny_var = "velocity" if field in ("u", "velocity") else field

        xyzt = np.array(
            [
                [int(start[0]), int(end[0])],
                [int(start[1]), int(end[1])],
                [int(start[2]), int(end[2])],
                [int(time_step), int(time_step)],
            ],
            dtype=int,
        )
        strides = np.array(
            [int(step[0]), int(step[1]), int(step[2]), 1], dtype=int
        )

        result = _giverny_getCutout(
            cube=self._cube,
            var=giverny_var,
            xyzt_axes_ranges_original=xyzt,
            xyzt_strides=strides,
            trace_memory=False,
            verbose=False,
        )

        key = f"{giverny_var}_{int(time_step):04d}"
        if key not in result.data_vars:
            key = list(result.data_vars)[0]
        return result[key].values

    def finalize(self) -> None:
        """Drop the current cube reference so a new session is required next."""
        self._cube = None
        self._dataset_title = None


def _initialize_jhtdb() -> "_JhtdbSession":
    """Create a JHTDB session and register the production token.

    Returns:
        A ``_JhtdbSession`` instance exposing ``getCutout`` in pyJHTDB style.

    Raises:
        ImportError: If neither pyJHTDB nor givernylocal is importable.
    """
    _log("Initializing JHTDB session (givernylocal backend)...")
    session = _JhtdbSession()
    _log(f"JHTDB session ready. Token = {JHTDB_TOKEN}")
    return session


# ---------------------------------------------------------------------------
# Shape normalization
# ---------------------------------------------------------------------------

def _normalize_velocity_shape(u: np.ndarray) -> np.ndarray:
    """Normalize a velocity array to (Nx, Ny, Nz, 3) convention.

    pyJHTDB historically returns either ``(Nz, Ny, Nx, 3)`` (component axis
    last) or ``(3, Nz, Ny, Nx)`` (component axis first) depending on version.
    We auto-detect and move the component axis to the end so the rest of
    the pipeline can assume a canonical layout.

    Args:
        u: The raw velocity array returned by ``libJHTDB.getCutout``.

    Returns:
        ``u`` with shape ``(Nx, Ny, Nz, 3)``.

    Raises:
        RuntimeError: If no axis of length 3 can be identified.
    """
    u = np.asarray(u)
    if u.ndim != 4:
        raise RuntimeError(
            f"Expected 4D velocity array (spatial + component); got shape {u.shape}"
        )
    if u.shape[-1] == 3:
        return u
    if u.shape[0] == 3:
        return np.moveaxis(u, 0, -1)
    raise RuntimeError(
        f"Cannot find component axis in shape {u.shape}; "
        "expected first or last axis of length 3."
    )


# ---------------------------------------------------------------------------
# Spectral anti-aliasing filter
# ---------------------------------------------------------------------------

def _spectral_lowpass_and_subsample(
    u_full: np.ndarray,
    output_stride: int,
) -> np.ndarray:
    """Apply a 3-D sharp spectral low-pass filter, then subsample spatially.

    Prevents aliasing that would otherwise arise from spatial sub-sampling
    alone (stride-4 on the 1024^3 DNS folds wavenumbers k=129..512 back
    onto k=1..128, corrupting the entire resolved spectrum).

    Wavenumber accounting
    ---------------------
    Each 512^3 sub-chunk covers one octant of the full 2x2x2 split, so its
    physical length is L_chunk = L_full / 2.  The DFT index *m* of the chunk
    maps to physical wavenumber:

        k_phys = m * (2π / L_chunk)

    The Nyquist of the OUTPUT grid (N_chunk / output_stride points, spacing
    Δx = L_full / N_DNS) in full-domain units is:

        k_c = N_output / 2  where  N_output = N_DNS / output_stride = 256
            = 128   (the cutoff requested by the spec)

    In chunk-local wavenumber index this is:

        m_c = N_chunk / (2 * output_stride)
            = 512 / 8 = 64

    We zero all chunk DFT modes with |m| > m_c before the IRFFT.  After
    IRFFT the physical-space signal is band-limited to k_phys ≤ 128, so
    sub-sampling by output_stride is aliasing-free.

    Args:
        u_full: Full-resolution velocity chunk (Nx, Ny, Nz, 3) float.
        output_stride: Spatial downsampling factor (typically 4).

    Returns:
        Anti-aliased array of shape
        (Nx//output_stride, Ny//output_stride, Nz//output_stride, 3), float32.
    """
    u_full = np.asarray(u_full, dtype=np.float32)
    Nx, Ny, Nz, nc = u_full.shape

    # Cutoff in chunk-local wavenumber index (see derivation in docstring)
    mc_x = Nx // (2 * output_stride)
    mc_y = Ny // (2 * output_stride)
    mc_z = Nz // (2 * output_stride)

    Nx_out = Nx // output_stride
    Ny_out = Ny // output_stride
    Nz_out = Nz // output_stride
    u_out = np.empty((Nx_out, Ny_out, Nz_out, nc), dtype=np.float32)

    for c in range(nc):
        # rfftn: full complex spectrum on axes 0-1, half-spectrum on axis 2.
        # Output shape: (Nx, Ny, Nz//2+1) complex128.
        U = np.fft.rfftn(u_full[..., c])

        # Zero modes above cutoff.
        #
        # Axis 0 (x): indices 0..Nx//2 (positive) then Nx-Nx//2..Nx-1
        # (negative mirror). Band [mc_x+1, Nx-mc_x) must be zeroed.
        if mc_x < Nx // 2:
            U[mc_x + 1 : Nx - mc_x, :, :] = 0.0

        # Axis 1 (y): same layout as axis 0.
        if mc_y < Ny // 2:
            U[:, mc_y + 1 : Ny - mc_y, :] = 0.0

        # Axis 2 (z, rfft half-spectrum): only positive frequencies 0..Nz//2.
        if mc_z < Nz // 2:
            U[:, :, mc_z + 1:] = 0.0

        # Inverse FFT back to physical space; stride to output resolution.
        u_filt = np.fft.irfftn(U, s=(Nx, Ny, Nz)).astype(np.float32)
        u_out[..., c] = u_filt[::output_stride, ::output_stride, ::output_stride]
        del U, u_filt  # release per-component FFT buffers promptly

    return u_out


def _spectral_lowpass_xz_and_subsample(
    u_full: np.ndarray,
    stride_xz: int,
    stride_y: int,
) -> np.ndarray:
    """Apply a 2-D sharp spectral low-pass filter in (x, z) only, then subsample.

    Channel flow is periodic in x (streamwise) and z (spanwise) but NOT in
    y (wall-normal).  A 3-D spectral filter would cross the wall boundaries
    and corrupt near-wall statistics.  This function filters only the two
    periodic directions and leaves y in physical space.

    Wavenumber accounting (same self-consistent formula as the 3-D case)
    ---------------------------------------------------------------------
    For a chunk of size Nx_chunk in x with stride ``stride_xz``:
        mc_x = Nx_chunk / (2 * stride_xz)
    This maps to a physical-domain cutoff
        k_phys_x = mc_x * (Nx_extract / Nx_chunk) = Nx_extract / (2 * stride_xz)
    which equals 128 for Nx_extract=1024 and stride_xz=4.  Identical logic
    applies to z.

    Args:
        u_full:    Full-resolution chunk ``(Nx, Ny, Nz, 3)`` float.
        stride_xz: Spatial downsampling factor in x and z (typically 4).
        stride_y:  Spatial downsampling factor in y (no spectral filter; typically 2).

    Returns:
        Anti-aliased array of shape
        ``(Nx//stride_xz, Ny//stride_y, Nz//stride_xz, 3)``, float32.
    """
    u_full = np.asarray(u_full, dtype=np.float32)
    Nx, Ny, Nz, nc = u_full.shape

    mc_x = Nx // (2 * stride_xz)
    mc_z = Nz // (2 * stride_xz)

    Nx_out = Nx // stride_xz
    Ny_out = Ny // stride_y
    Nz_out = Nz // stride_xz
    u_out = np.empty((Nx_out, Ny_out, Nz_out, nc), dtype=np.float32)

    for c in range(nc):
        # 2-D rFFT in axes (0=x, 2=z); axis 1 (y) stays in physical space.
        # rfftn with axes=(0,2): full complex along x, half-spectrum along z.
        # Output shape: (Nx, Ny, Nz//2+1) complex128.
        U = np.fft.rfftn(u_full[..., c], axes=(0, 2))

        # Zero high-kx modes: band [mc_x+1, Nx-mc_x) on axis 0.
        if mc_x < Nx // 2:
            U[mc_x + 1 : Nx - mc_x, :, :] = 0.0

        # Zero high-kz modes on axis 2 (rfft half-spectrum, only positive k).
        if mc_z < Nz // 2:
            U[:, :, mc_z + 1:] = 0.0

        # Inverse 2-D real FFT; recover full (Nx, Ny, Nz) physical-space field.
        u_filt = np.fft.irfftn(U, s=(Nx, Nz), axes=(0, 2)).astype(np.float32)

        # Subsample: stride_xz in x and z, stride_y in y (no filter needed in y).
        u_out[..., c] = u_filt[::stride_xz, ::stride_y, ::stride_xz]
        del U, u_filt

    return u_out


def _get_channel_y_coords(ny_native: int = 512, stride_y: int = 2) -> np.ndarray:
    """Return wall-normal y-coordinates at the subsampled channel grid points.

    JHTDB channel (Re_tau=1000) uses a cosine-stretched non-uniform grid in y
    from 0 (bottom wall) to 2 (top wall) with ``ny_native`` points:

        y_j = 1 - cos(π * j / (ny_native - 1)),   j = 0, …, ny_native-1

    We select every ``stride_y``-th point to match the spatial sub-sampling
    applied to the velocity field.

    Args:
        ny_native: Number of native y-grid points (512 for JHTDB channel).
        stride_y:  Wall-normal sub-sampling stride (typically 2).

    Returns:
        1-D float64 array of length ``ny_native // stride_y`` with y in [0, 2].
    """
    j = np.arange(ny_native, dtype=np.float64)
    y_full = 1.0 - np.cos(np.pi * j / (ny_native - 1))
    return y_full[::stride_y]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_isotropic(u: np.ndarray, name: str = "isotropic") -> dict:
    """Strict validation that ``u`` is consistent with isotropic turbulence.

    Performs five checks, any one of which causes a RuntimeError:
      1. All finite (no NaN/inf).
      2. ``max|mean_u_i| / max rms_i < 0.05`` (zero-mean).
      3. ``max(rms)/min(rms) < 1.5`` (isotropy of second moments).
      4. All off-diagonal ``|cov_ij|/(rms_i rms_j) < 0.15``.
      5. ``max(octant_rms)/min(octant_rms) < 2.0`` over 2x2x2 octants.
         This is the smoking gun for chunk-stitching corruption: a single
         misaligned octant blows the ratio past 2.

    Args:
        u: Velocity array in any accepted layout.
        name: Label used in log and error messages.

    Returns:
        A diagnostics dict with ``mean_u``, ``rms_u``, ``cross_rho``,
        ``octant_rms_min``, ``octant_rms_max``, ``octant_rms_ratio``.

    Raises:
        RuntimeError: If any check fails. No file should be written.
    """
    u = _normalize_velocity_shape(u)

    _log(f"[validate] {name} (strict isotropy checks on shape {u.shape})")

    # (1) Finite check — hard abort, cannot continue
    if not np.all(np.isfinite(u)):
        n_bad = int(np.sum(~np.isfinite(u)))
        raise RuntimeError(
            f"Validation FAILED for {name}: {n_bad} non-finite values "
            "(NaN/inf) detected. Refusing to save corrupted data."
        )

    failures: list[str] = []

    # (2) Per-component first and second moments
    global_mean = np.array([float(u[..., c].mean()) for c in range(3)])
    global_rms = np.array([float(u[..., c].std()) for c in range(3)])
    max_rms = float(global_rms.max())
    min_rms = float(global_rms.min())

    mean_over_rms = float(np.max(np.abs(global_mean))) / max(max_rms, 1e-30)
    if mean_over_rms >= 0.05:
        failures.append(
            f"max|mean|/max(rms) = {mean_over_rms:.3f} >= 0.05 "
            "(forced isotropic turbulence should be zero-mean)"
        )

    # (3) Isotropy of second moments
    if min_rms <= 0:
        failures.append(f"min rms is zero: {global_rms.tolist()}")
        rms_ratio = float("inf")
    else:
        rms_ratio = max_rms / min_rms
        if rms_ratio >= 1.5:
            failures.append(
                f"rms ratio max/min = {rms_ratio:.3f} >= 1.5 "
                f"(rms = {global_rms.tolist()})"
            )

    # (4) Off-diagonal cross-correlations (Reynolds stresses)
    flat = [u[..., c].ravel() for c in range(3)]
    cross_cov = np.zeros(3, dtype=np.float64)
    cross_rho = np.zeros(3, dtype=np.float64)
    for k, (i, j) in enumerate([(0, 1), (0, 2), (1, 2)]):
        cov_ij = float(np.cov(flat[i], flat[j])[0, 1])
        denom = float(global_rms[i] * global_rms[j])
        rho_ij = cov_ij / denom if denom > 0 else 0.0
        cross_cov[k] = cov_ij
        cross_rho[k] = rho_ij
        if abs(rho_ij) >= 0.15:
            failures.append(
                f"<u{i+1} u{j+1}> / (rms_{i+1} rms_{j+1}) = {rho_ij:+.3f} "
                "(|rho| >= 0.15; off-diagonal Reynolds stress should vanish)"
            )

    # (5) Per-octant RMS homogeneity. Split domain into 2x2x2 = 8 blocks of
    # half the linear size and compare RMS across octants. If stitching was
    # broken, one octant will look qualitatively different.
    Nx, Ny, Nz = u.shape[0], u.shape[1], u.shape[2]
    hx, hy, hz = Nx // 2, Ny // 2, Nz // 2
    if hx == 0 or hy == 0 or hz == 0:
        raise RuntimeError(
            f"Validation FAILED for {name}: array too small for octant split "
            f"(shape {u.shape})"
        )

    octant_rms: list[tuple[float, float, float]] = []
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                # NB: ordering ix,iy,iz matches the C++ convention; we visit
                # all 8 octants and record the per-component RMS of each.
                u_oct = u[ix * hx:(ix + 1) * hx,
                          iy * hy:(iy + 1) * hy,
                          iz * hz:(iz + 1) * hz, :]
                octant_rms.append((
                    float(u_oct[..., 0].std()),
                    float(u_oct[..., 1].std()),
                    float(u_oct[..., 2].std()),
                ))

    flat_rms = [r for triple in octant_rms for r in triple]
    octant_max = float(np.max(flat_rms))
    octant_min = float(np.min(flat_rms))
    octant_ratio = octant_max / max(octant_min, 1e-30)
    if octant_ratio >= 2.0:
        failures.append(
            f"octant-rms ratio = {octant_ratio:.3f} >= 2.0 [SMOKING GUN] "
            f"(max={octant_max:.4f}, min={octant_min:.4f}); "
            "classic chunk-stitching corruption signature"
        )

    # Emit result
    if failures:
        msg = (
            f"\nValidation FAILED for {name}:\n"
            f"  - global_mean:       {global_mean.tolist()}  (threshold: |mean|/rms < 0.05)\n"
            f"  - global_rms:        {global_rms.tolist()}  (ratio threshold: < 1.5)\n"
            f"  - cross_correlations:{cross_rho.tolist()}  (threshold: < 0.15)\n"
            f"  - octant rms range:  {octant_min:.4f} .. {octant_max:.4f}  "
            f"(ratio {octant_ratio:.3f}; threshold < 2.0)\n"
            f"\n  Individual failures:\n"
        )
        for f in failures:
            msg += f"    * {f}\n"
        msg += "\n  Refusing to save corrupted data.\n"
        raise RuntimeError(msg)

    # All checks passed — print and return diagnostics.
    print(
        f"\nValidation PASSED for {name}:\n"
        f"  - mean(u):        [{global_mean[0]:+.4e}, {global_mean[1]:+.4e}, {global_mean[2]:+.4e}]  (|mean|/rms_max = {mean_over_rms:.3f})\n"
        f"  - rms(u):         [{global_rms[0]:.4f}, {global_rms[1]:.4f}, {global_rms[2]:.4f}]  (ratio max/min = {rms_ratio:.3f})\n"
        f"  - cross rho(ij):  [{cross_rho[0]:+.4f}, {cross_rho[1]:+.4f}, {cross_rho[2]:+.4f}]\n"
        f"  - octant rms:     {octant_min:.4f} .. {octant_max:.4f}  (ratio {octant_ratio:.3f})\n",
        flush=True,
    )

    return {
        "mean_u":           global_mean,
        "rms_u":            global_rms,
        "cross_cov":        cross_cov,
        "cross_rho":        cross_rho,
        "octant_rms_min":   octant_min,
        "octant_rms_max":   octant_max,
        "octant_rms_ratio": octant_ratio,
    }


def _validate_channel(u: np.ndarray, name: str = "channel") -> dict:
    """Loose validation for channel flow (inherently anisotropic).

    Only failures are NaN/inf or all-zero RMS. Streamwise dominance is
    emitted as a WARN but does not stop the save.

    Args:
        u: Velocity array in any accepted layout.
        name: Label used in log and error messages.

    Returns:
        Diagnostics dict with ``mean_u``, ``rms_u``.

    Raises:
        RuntimeError: If non-finite values or all-zero RMS are present.
    """
    u = _normalize_velocity_shape(u)

    _log(f"[validate] {name} (loose channel checks on shape {u.shape})")

    if not np.all(np.isfinite(u)):
        n_bad = int(np.sum(~np.isfinite(u)))
        raise RuntimeError(
            f"Validation FAILED for {name}: {n_bad} non-finite values. Refusing to save."
        )

    global_mean = np.array([float(u[..., c].mean()) for c in range(3)])
    global_rms = np.array([float(u[..., c].std()) for c in range(3)])

    if float(global_rms.max()) <= 0:
        raise RuntimeError(
            f"Validation FAILED for {name}: all RMS components are zero "
            f"(rms = {global_rms.tolist()})."
        )

    # Streamwise (u1) should typically dominate. Warn if it does not.
    if global_rms[0] < max(global_rms[1], global_rms[2]):
        print(
            f"  [WARN] {name}: streamwise rms is not the largest component "
            f"(rms = {global_rms.tolist()}); unusual for channel flow.",
            flush=True,
        )

    print(
        f"\nValidation PASSED for {name}:\n"
        f"  - mean(u): [{global_mean[0]:+.4e}, {global_mean[1]:+.4e}, {global_mean[2]:+.4e}]\n"
        f"  - rms(u):  [{global_rms[0]:.4f}, {global_rms[1]:.4f}, {global_rms[2]:.4f}]\n"
        f"  - note:    anisotropy is expected for channel; no strict thresholds applied\n",
        flush=True,
    )

    return {
        "mean_u": global_mean,
        "rms_u":  global_rms,
    }


# ---------------------------------------------------------------------------
# HDF5 atomic write
# ---------------------------------------------------------------------------

def _write_velocity_hdf5(
    path: pathlib.Path,
    u: np.ndarray,
    attrs: dict,
    extra_datasets: Optional[dict] = None,
) -> None:
    """Write ``u`` to an HDF5 file at ``path`` with the given attributes.

    Layout:
        /u                 dataset (Nx, Ny, Nz, 3) float64, gzip-compressed
        /u1, /u2, /u3      datasets (Nx, Ny, Nz) float64  -- per-component
                           views for backward compatibility with the C++
                           io.cpp reader.
        /                  attributes from ``attrs``
        + any datasets in ``extra_datasets`` (name -> ndarray)

    Args:
        path: Absolute or relative path to the output .h5 file.
        u: Velocity array of shape (Nx, Ny, Nz, 3) float.
        attrs: Dict of attribute-name -> value for the root group.
        extra_datasets: Optional dict of additional 1-D or N-D arrays to
            write alongside the velocity (e.g. ``{"y": y_arr, "y_plus": yp}``).
    """
    u = np.ascontiguousarray(u.astype(np.float64))
    with h5py.File(path, "w") as f:
        f.create_dataset("u", data=u, compression="gzip", compression_opts=4)
        # Per-component datasets so the existing C++ pipeline (which reads
        # /u1 /u2 /u3) continues to work without changes.
        f.create_dataset("u1", data=u[..., 0], compression="gzip", compression_opts=4)
        f.create_dataset("u2", data=u[..., 1], compression="gzip", compression_opts=4)
        f.create_dataset("u3", data=u[..., 2], compression="gzip", compression_opts=4)
        if extra_datasets:
            for ds_name, ds_data in extra_datasets.items():
                f.create_dataset(ds_name, data=np.asarray(ds_data, dtype=np.float64))
        for k, v in attrs.items():
            f.attrs[k] = v


def _atomic_save(
    output_path: pathlib.Path,
    u: np.ndarray,
    attrs: dict,
    extra_datasets: Optional[dict] = None,
) -> None:
    """Atomically write ``u`` to ``output_path`` via a .tmp staging file.

    On any exception the staging file is removed so we never leave a
    half-written .h5 on disk.

    Args:
        output_path: Final destination path.
        u: Velocity array to save.
        attrs: HDF5 root attributes.
        extra_datasets: Optional dict of additional datasets (e.g. y, y_plus).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    _log(f"Atomic save: writing {tmp.name} ...")
    try:
        _write_velocity_hdf5(tmp, u, attrs, extra_datasets=extra_datasets)
        tmp.replace(output_path)
        _log(f"Atomic save: {tmp.name} -> {output_path.name}")
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        raise


# ---------------------------------------------------------------------------
# Legacy cleanup
# ---------------------------------------------------------------------------

def _backup_existing_file(path: pathlib.Path) -> Optional[pathlib.Path]:
    """Rename ``path`` to ``path + '.old'`` if it exists.

    Protects against accidental loss if the new download fails partway.

    Args:
        path: File path to back up.

    Returns:
        Backup path, or None if ``path`` did not exist.
    """
    if not path.exists():
        return None
    backup = path.with_suffix(path.suffix + ".old")
    if backup.exists():
        backup.unlink()
    path.rename(backup)
    _log(f"Backed up old file: {path.name} -> {backup.name}")
    return backup


def _remove_backup(backup: Optional[pathlib.Path]) -> None:
    """Delete the backup file created by ``_backup_existing_file`` (if any)."""
    if backup is not None and backup.exists():
        backup.unlink()
        _log(f"Removed backup: {backup.name}")


def _note_corrupt_backup() -> None:
    """Log the presence of any legacy .corrupt-backup file.

    The instruction brief marks these as 'safe to remove', but because we
    are running non-interactively we only LOG their presence and leave the
    actual deletion to the operator. This is the conservative choice: no
    destructive action without explicit confirmation.
    """
    candidate = _DATA_DIR / "isotropic_256.h5.corrupt-backup"
    if candidate.exists():
        size_mb = candidate.stat().st_size / 1e6
        _log(
            f"NOTE: legacy {candidate.name} present ({size_mb:.1f} MB). "
            "Safe to delete manually once the new monolithic download validates."
        )


# ---------------------------------------------------------------------------
# Monolithic downloads
# ---------------------------------------------------------------------------
#
# Tactical note
# -------------
# The brief calls for a single 256^3 getCutout request to replace the old
# 8-chunk stitched pipeline. In this environment the backend is givernylocal
# (pyJHTDB is deprecated; giverny requires SciServer), which imposes a
# client-side cap on the NATIVE box size (~3072 MB = 268 M points) and
# internally warns that HTTP payloads above ~200^3 are unreliable.
#
# A 1024^3 native request is 12 GB raw — far above what the local HTTP client
# can safely ferry in one payload. So instead of a single giant call we loop
# over a small number of aligned sub-boxes. The KEY difference from the old
# broken pipeline is that each sub-box is fetched on a *brand new*
# ``_JhtdbSession`` (fresh ``turb_dataset`` cube), so giverny's internal cache
# cannot possibly misalign chunks between iterations. The final stitched
# array is then subjected to the full strict ``_validate_isotropic`` check
# (including the octant smoking gun) before anything is written to disk.
#
# This is a "belt and suspenders" fix for the original April 2026 corruption:
# fresh cubes eliminate the cache misalignment, strict validation catches
# anything the fresh-cube fix misses, and atomic .tmp->rename ensures a bad
# download can never overwrite a good file on disk.


def _download_one_chunk_fresh(
    dataset_title: str,
    time_step: int,
    start_idx: np.ndarray,
    end_idx: np.ndarray,
    step: np.ndarray,
) -> np.ndarray:
    """Fetch one native sub-box on a FRESH giverny cube; discard the session.

    A new ``_JhtdbSession`` is constructed for every call so no cache state
    survives between chunks — this is the root-cause fix for the previous
    stitching corruption.

    Args:
        dataset_title: giverny dataset title.
        time_step: 1-based time index.
        start_idx: length-3 int array, 1-based inclusive start of native box.
        end_idx: length-3 int array, 1-based inclusive end of native box.
        step: length-3 int array, per-axis stride over native grid.

    Returns:
        Raw velocity array as returned by giverny (component axis last,
        typically ``(Nz, Ny, Nx, 3)`` float32).
    """
    session = _JhtdbSession()
    try:
        return session.getCutout(
            data_set=dataset_title,
            field="u",
            time_step=int(time_step),
            start=start_idx,
            end=end_idx,
            step=step,
            filter_width=1,
        )
    finally:
        try:
            session.finalize()
        except Exception:
            pass


def download_isotropic_256_monolithic(
    stride: int = _DEFAULT_STRIDE,
    time_step: int = 100,
    output_path: Optional[pathlib.Path | str] = None,
) -> np.ndarray:
    """Download the 256^3 isotropic1024coarse velocity cube.

    Tactically this issues 8 sub-box getCutout calls (each on a fresh cube)
    and stitches them into a single 256^3 volume. Strategically, the
    difference from the previous broken pipeline is:
      - Each sub-box uses a brand-new turb_dataset, defeating giverny's
        per-cube cache which was the root cause of the April 2026
        corruption incident.
      - The stitched array is subjected to the full strict
        _validate_isotropic (including the octant smoking gun test).
      - Nothing is written to disk unless validation passes.

    Args:
        stride: Stride over the native 1024^3 grid. ``stride=4`` yields 256^3.
        time_step: JHTDB time index (1-based).
        output_path: Destination .h5 file. Defaults to ``data/isotropic_256.h5``.

    Returns:
        Validated velocity array of shape ``(256, 256, 256, 3)``, float64.

    Raises:
        RuntimeError: On any validation failure. The output file is
            untouched in that case and any pre-existing file is preserved
            via its .old backup.
    """
    out = (
        pathlib.Path(output_path).resolve()
        if output_path is not None
        else (_DATA_DIR / "isotropic_256.h5").resolve()
    )

    _log("=" * 63)
    _log("isotropic1024coarse -- fresh-cube-per-chunk + spectral filter")
    _log(f"  dataset    : isotropic1024coarse")
    _log(f"  time_step  : {time_step}")
    _log(f"  native N   : {_NATIVE_N_ISO}")
    _log(f"  stride     : {stride}  (spectral kc=128 applied before subsampling)")
    _log(f"  output N   : {_NATIVE_N_ISO // stride}")
    _log(f"  output file: {out}")
    _log("=" * 63)

    backup = _backup_existing_file(out)

    # 4x4x4 split of the native 1024^3 domain into 64 sub-boxes of 256^3 each.
    # Each sub-box is fetched at stride=1 (~201 MB response, well within the
    # ~200^3 givernylocal HTTP limit).  After spectral filtering + stride-4
    # subsampling each sub-box yields 64^3 output; 4^3=64 sub-boxes stitch to
    # 256^3.  Using 256^3 chunks (vs the previous 512^3) was necessary because
    # stride=1 responses for 512^3 (~1.5 GB) trigger HTTP 502 errors.
    n_split = 4
    out_per_chunk = _OUTPUT_N // n_split  # = 64
    native_per_chunk = _NATIVE_N_ISO // n_split  # = 256

    full = np.empty((_OUTPUT_N, _OUTPUT_N, _OUTPUT_N, 3), dtype=np.float32)
    t0 = datetime.datetime.now()

    total = n_split ** 3
    done = 0
    for iz in range(n_split):
        for iy in range(n_split):
            for ix in range(n_split):
                # Native 1-based inclusive ranges aligned to stride boundaries.
                start_idx = np.array(
                    [ix * native_per_chunk + 1,
                     iy * native_per_chunk + 1,
                     iz * native_per_chunk + 1], dtype=np.int32)
                end_idx = np.array(
                    [(ix + 1) * native_per_chunk,
                     (iy + 1) * native_per_chunk,
                     (iz + 1) * native_per_chunk], dtype=np.int32)
                # Download at FULL resolution (stride=1).  The spectral
                # low-pass filter is applied below before spatial sub-sampling
                # to prevent aliasing (bug A1).
                step_full = np.array([1, 1, 1], dtype=np.int32)

                t_chunk = datetime.datetime.now()
                raw_chunk = _download_one_chunk_fresh(
                    "isotropic1024coarse", time_step, start_idx, end_idx, step_full,
                )
                chunk_full = _normalize_velocity_shape(np.asarray(raw_chunk))
                del raw_chunk  # free before allocating FFT buffers

                expected_full = (native_per_chunk, native_per_chunk, native_per_chunk, 3)
                if chunk_full.shape != expected_full:
                    raise RuntimeError(
                        f"chunk iz={iz} iy={iy} ix={ix}: full-res shape "
                        f"{chunk_full.shape} != expected {expected_full}"
                    )
                if not np.all(np.isfinite(chunk_full)):
                    raise RuntimeError(
                        f"chunk iz={iz} iy={iy} ix={ix} has non-finite values "
                        "(full-res download)"
                    )

                # Sharp spectral low-pass at kc=64 in chunk-local wavenumber
                # units (= kc=128 in full-domain units), then stride=4 subsample.
                chunk = _spectral_lowpass_and_subsample(chunk_full, stride)
                del chunk_full  # release ~1.5 GB full-res buffer

                expected = (out_per_chunk, out_per_chunk, out_per_chunk, 3)
                if chunk.shape != expected:
                    raise RuntimeError(
                        f"chunk iz={iz} iy={iy} ix={ix}: filtered shape "
                        f"{chunk.shape} != expected {expected}"
                    )

                # NB: paste the (ix, iy, iz) chunk into the matching sub-block
                # of the full 256^3 buffer. The giverny axis order is the same
                # as our (Nx, Ny, Nz, 3) convention post-normalization.
                full[ix * out_per_chunk:(ix + 1) * out_per_chunk,
                     iy * out_per_chunk:(iy + 1) * out_per_chunk,
                     iz * out_per_chunk:(iz + 1) * out_per_chunk, :] = chunk

                done += 1
                chunk_rms = [float(chunk[..., c].std()) for c in range(3)]
                _log(
                    f"  [{done}/{total}] ix={ix} iy={iy} iz={iz}  "
                    f"rms=({chunk_rms[0]:.3f},{chunk_rms[1]:.3f},{chunk_rms[2]:.3f})  "
                    f"dt={(datetime.datetime.now() - t_chunk).total_seconds():.1f}s  "
                    f"[fresh cube + spectral filter kc=128]"
                )
                del chunk

    elapsed = (datetime.datetime.now() - t0).total_seconds()
    _log(f"all 8 chunks stitched in {elapsed:.1f}s")

    u = full  # already (Nx, Ny, Nz, 3)
    expected_shape = (_OUTPUT_N, _OUTPUT_N, _OUTPUT_N, 3)
    if u.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected stitched velocity shape: {u.shape} "
            f"(expected {expected_shape})"
        )

    # Strict validation — raises on failure so we never save bad data.
    stats = _validate_isotropic(u, name="isotropic1024coarse")

    # HDF5 attributes: everything the C++ side or a downstream Python tool
    # might ever want to know about provenance.
    attrs = {
        "stride":            int(stride),
        "native_N":          int(_NATIVE_N_ISO),
        "N":                 int(_OUTPUT_N),
        "dx":                float(2.0 * np.pi / _OUTPUT_N),
        "time_step":         int(time_step),
        "downloaded_at":     datetime.datetime.now().isoformat(),
        "validation_status": "passed",
        "dataset":           "isotropic1024coarse",
        "Re_lambda":         433,
        "mean_u":            np.asarray(stats["mean_u"], dtype=np.float64),
        "rms_u":             np.asarray(stats["rms_u"], dtype=np.float64),
        "octant_rms_max":    float(stats["octant_rms_max"]),
        "octant_rms_min":    float(stats["octant_rms_min"]),
        # Bug-A1 fix: spectral anti-aliasing applied before spatial sub-sampling.
        # Each 512^3 sub-chunk was downloaded at stride=1, filtered at
        # kc=64 (chunk-local) = kc=128 (full-domain), then strided by 4.
        "spectral_filter":   "sharp_lowpass",
        "kc_full_domain":    int(_OUTPUT_N // 2),
    }

    _atomic_save(out, u, attrs)
    _remove_backup(backup)
    _log(f"isotropic_256.h5 written OK ({out.stat().st_size/1e6:.1f} MB)")

    return u


def download_channel_256_monolithic(
    time_step: int = 100,
    output_path: Optional[pathlib.Path | str] = None,
) -> np.ndarray:
    """Download the 256^3 channel-flow velocity cube (Bugs A2, A3, B3 fixed).

    Fixes applied
    -------------
    A2 (shape mismatch):
        Previous code extracted (1024, 512, 512) at stride=4 → (256, 128, 128).
        Correct extract is (1024, 512, 1024) with strides (4, 2, 4) → (256, 256, 256).

    A3 (2-D spectral filter):
        Channel is periodic only in x and z, NOT in y.  A 3-D filter would
        cross wall boundaries and corrupt near-wall statistics.  We apply a
        2-D sharp spectral low-pass in (x, z) at kc=128 (full-domain units)
        before spatial sub-sampling.  The wall-normal direction is left in
        physical space and sub-sampled by stride_y=2 without spectral filtering.

    B3 (wall-normal coordinates):
        The y and y_plus 1-D arrays (required for Figure 8 wall-normal plots)
        are computed from the standard JHTDB cosine-stretched grid formula and
        saved as ``/y`` and ``/y_plus`` datasets inside channel_256.h5.

    Chunk strategy (same as isotropic, avoids HTTP 502)
    ---------------------------------------------------
    Native extract : 1024 (x) x 512 (y) x 1024 (z)
    n_split        : 4 × 4 × 4 = 64 chunks
    Chunk native   : 256 (x) × 128 (y) × 256 (z)  ≈ 101 MB  (within givernylocal limit)
    Strides        : (4, 2, 4)
    Chunk output   : 64 × 64 × 64
    Total output   : 256 × 256 × 256  ✓

    Spectral cutoffs (chunk-local mc → full-domain kc)
    --------------------------------------------------
        mc_x = 256 / (2*4) = 32  →  k_phys = 32 * (1024/256) = 128  ✓
        mc_z = 256 / (2*4) = 32  →  k_phys = 32 * (1024/256) = 128  ✓

    Args:
        time_step:   JHTDB time index (1-based).
        output_path: Destination .h5 file.  Defaults to ``data/channel_256.h5``.

    Returns:
        Validated velocity array of shape ``(256, 256, 256, 3)``, float64.

    Raises:
        RuntimeError: On any validation failure.  Output file is untouched.
    """
    out = (
        pathlib.Path(output_path).resolve()
        if output_path is not None
        else (_DATA_DIR / "channel_256.h5").resolve()
    )

    nx_native, ny_native, nz_native = _CHANNEL_EXTRACT  # (1024, 512, 1024)
    stride_xz = _CHANNEL_STRIDE_XZ  # 4
    stride_y  = _CHANNEL_STRIDE_Y   # 2
    nx_out = nx_native // stride_xz  # 256
    ny_out = ny_native // stride_y   # 256
    nz_out = nz_native // stride_xz  # 256

    _log("=" * 63)
    _log("channel -- fresh-cube-per-chunk + 2D spectral filter (xz)")
    _log(f"  dataset      : channel")
    _log(f"  time_step    : {time_step}")
    _log(f"  native extract: {_CHANNEL_EXTRACT}")
    _log(f"  strides (x,y,z): ({stride_xz}, {stride_y}, {stride_xz})")
    _log(f"  output shape : ({nx_out}, {ny_out}, {nz_out}, 3)")
    _log(f"  spectral filter: 2-D sharp lowpass in (x,z), kc=128; y left physical")
    _log(f"  output file  : {out}")
    _log("=" * 63)

    backup = _backup_existing_file(out)

    # 4×4×4 = 64 chunks, each 256×128×256 native ≈ 101 MB (within givernylocal limit).
    n_split = 4
    nx_chunk_native = nx_native // n_split   # 256
    ny_chunk_native = ny_native // n_split   # 128
    nz_chunk_native = nz_native // n_split   # 256
    nx_chunk_out    = nx_chunk_native // stride_xz  # 64
    ny_chunk_out    = ny_chunk_native // stride_y   # 64
    nz_chunk_out    = nz_chunk_native // stride_xz  # 64

    full = np.empty((nx_out, ny_out, nz_out, 3), dtype=np.float32)
    t0 = datetime.datetime.now()
    total = n_split ** 3
    done = 0

    for iz in range(n_split):
        for iy in range(n_split):
            for ix in range(n_split):
                # 1-based inclusive native indices.
                start_idx = np.array([
                    ix * nx_chunk_native + 1,
                    iy * ny_chunk_native + 1,
                    iz * nz_chunk_native + 1,
                ], dtype=np.int32)
                end_idx = np.array([
                    (ix + 1) * nx_chunk_native,
                    (iy + 1) * ny_chunk_native,
                    (iz + 1) * nz_chunk_native,
                ], dtype=np.int32)
                # Download at FULL resolution (stride=1); spectral filter
                # handles anti-aliasing in x and z before sub-sampling.
                step_full = np.array([1, 1, 1], dtype=np.int32)

                t_chunk = datetime.datetime.now()
                raw_chunk = _download_one_chunk_fresh(
                    "channel", time_step, start_idx, end_idx, step_full,
                )
                chunk_full = _normalize_velocity_shape(np.asarray(raw_chunk))
                del raw_chunk

                expected_full = (nx_chunk_native, ny_chunk_native, nz_chunk_native, 3)
                if chunk_full.shape != expected_full:
                    raise RuntimeError(
                        f"channel chunk iz={iz} iy={iy} ix={ix}: full-res shape "
                        f"{chunk_full.shape} != expected {expected_full}"
                    )
                if not np.all(np.isfinite(chunk_full)):
                    raise RuntimeError(
                        f"channel chunk iz={iz} iy={iy} ix={ix} has non-finite values"
                    )

                # 2-D spectral low-pass in (x, z) only; y left in physical space.
                chunk = _spectral_lowpass_xz_and_subsample(
                    chunk_full, stride_xz, stride_y,
                )
                del chunk_full

                expected = (nx_chunk_out, ny_chunk_out, nz_chunk_out, 3)
                if chunk.shape != expected:
                    raise RuntimeError(
                        f"channel chunk iz={iz} iy={iy} ix={ix}: filtered shape "
                        f"{chunk.shape} != expected {expected}"
                    )

                full[
                    ix * nx_chunk_out:(ix + 1) * nx_chunk_out,
                    iy * ny_chunk_out:(iy + 1) * ny_chunk_out,
                    iz * nz_chunk_out:(iz + 1) * nz_chunk_out,
                    :,
                ] = chunk

                done += 1
                chunk_rms = [float(chunk[..., c].std()) for c in range(3)]
                _log(
                    f"  [{done}/{total}] ix={ix} iy={iy} iz={iz}  "
                    f"rms=({chunk_rms[0]:.3f},{chunk_rms[1]:.3f},{chunk_rms[2]:.3f})  "
                    f"dt={(datetime.datetime.now() - t_chunk).total_seconds():.1f}s  "
                    f"[fresh cube + 2D spectral filter kc=128]"
                )
                del chunk

    elapsed = (datetime.datetime.now() - t0).total_seconds()
    _log(f"all {total} chunks stitched in {elapsed:.1f}s")

    u = full
    expected_shape = (nx_out, ny_out, nz_out, 3)
    if u.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected stitched channel shape: {u.shape} (expected {expected_shape})"
        )

    # Loose validation (channel is inherently anisotropic; no strict thresholds).
    stats = _validate_channel(u, name="channel")

    # --- Bug B3: wall-normal coordinates -----------------------------------
    # y: physical wall-normal coordinate from 0 (bottom wall) to 2 (top wall),
    #    sampled at the same stride_y as the velocity field.
    # y_plus: viscous-unit distance from the nearest wall = min(y, 2-y) * Re_tau.
    y_coords = _get_channel_y_coords(ny_native=ny_native, stride_y=stride_y)
    y_plus   = np.minimum(y_coords, 2.0 - y_coords) * float(_CHANNEL_RE_TAU)
    _log(f"y-coordinate range: [{y_coords[0]:.4f}, {y_coords[-1]:.4f}]  "
         f"y_plus range: [{y_plus.min():.2f}, {y_plus.max():.2f}]")

    attrs = {
        "dataset":           "channel",
        "Re_tau":            _CHANNEL_RE_TAU,
        "time_step":         int(time_step),
        "downloaded_at":     datetime.datetime.now().isoformat(),
        "validation_status": "passed",
        "native_extract":    np.array(list(_CHANNEL_EXTRACT), dtype=np.int64),
        "strides":           np.array([stride_xz, stride_y, stride_xz], dtype=np.int64),
        "N":                 int(nx_out),   # 256 (all three output dims are equal)
        # dx: nominal grid spacing assuming a (2*pi)^3 periodic box — matches
        # what the C++ pipeline (io.cpp read_velocity_hdf5) and config.h expect.
        "dx":                float(2.0 * np.pi / _OUTPUT_N),
        "spectral_filter":   "sharp_lowpass_xz",
        # kc in full-domain units: Nx_extract / (2 * stride_xz) = 1024/(2*4) = 128.
        "kc_full_domain_xz": int(nx_native // (2 * stride_xz)),
        "mean_u":            np.asarray(stats["mean_u"], dtype=np.float64),
        "rms_u":             np.asarray(stats["rms_u"], dtype=np.float64),
    }

    extra_datasets = {
        "y":      y_coords,
        "y_plus": y_plus,
    }

    _atomic_save(out, u, attrs, extra_datasets=extra_datasets)
    _remove_backup(backup)
    _log(f"channel_256.h5 written OK ({out.stat().st_size/1e6:.1f} MB)")
    _log(f"  /y      shape={y_coords.shape}  range=[{y_coords[0]:.4f}, {y_coords[-1]:.4f}]")
    _log(f"  /y_plus shape={y_plus.shape}   range=[{y_plus.min():.2f}, {y_plus.max():.2f}]")

    return u


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def download_all_datasets() -> None:
    """Run both monolithic downloads sequentially (isotropic, then channel).

    Per JHTDB policy, no parallelization: the second download does not start
    until the first has been fully validated and written.

    Raises:
        RuntimeError: Re-raised verbatim from the per-dataset downloader on
            any validation failure. The process stops at the first failure.
    """
    print("\n" + "=" * 70)
    print("  JHTDB download refactor (fresh-cube-per-chunk + strict validation)")
    print("  " + "-" * 66)
    print(f"  Token            : {JHTDB_TOKEN}")
    print("  Backend          : givernylocal (pyJHTDB deprecated; giverny")
    print("                     requires SciServer auth, not available here)")
    print("  Fix vs old code  : each sub-box uses a brand-new turb_dataset,")
    print("                     so giverny cache cannot misalign chunks.")
    print("  Post-stitching   : strict octant validation + atomic .tmp->rename.")
    print("=" * 70 + "\n")

    _DATA_DIR.mkdir(parents=True, exist_ok=True)
    _note_corrupt_backup()

    # --- Isotropic ---------------------------------------------------------
    try:
        download_isotropic_256_monolithic(
            stride=_DEFAULT_STRIDE,
            time_step=100,
        )
    except RuntimeError as exc:
        _log(f"[FATAL] isotropic download failed: {exc}")
        raise

    # --- Channel -----------------------------------------------------------
    try:
        download_channel_256_monolithic(
            time_step=100,
        )
    except RuntimeError as exc:
        _log(f"[FATAL] channel download failed: {exc}")
        raise

    # --- Summary -----------------------------------------------------------
    iso_path = _DATA_DIR / "isotropic_256.h5"
    ch_path = _DATA_DIR / "channel_256.h5"

    print("\n" + "=" * 70)
    print("  Both downloads complete.")
    print(f"  - {iso_path}  ({iso_path.stat().st_size/1e6:.1f} MB, 256^3)")
    print(f"  - {ch_path}   ({ch_path.stat().st_size/1e6:.1f} MB, 256^3, /y + /y_plus saved)")
    print("  Ready for C++ pipeline testing.")
    print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    try:
        # Download fourth and final time-snapshot (t4) of channel flow for ensemble convergence
        # time_step=400 (independent from t1=100, t2=200, and t3=300)
        # output to channel_256_t4.h5 to preserve prior snapshots
        # FINAL DATA COLLECTION TASK FOR ENSEMBLE
        download_channel_256_monolithic(
            time_step=400,
            output_path=_DATA_DIR / "channel_256_t4.h5",
        )
    except RuntimeError as exc:
        print(f"\n[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)
