"""
FractalToroid — continuum-field simulation of Bob Greenyer's "fractal
toroidal moment" (see README.md). Dimensionless plasma units throughout:
c = ε0 = μ0 = 1, lengths in electron skin depths, times in 1/ω_p.

Phase 1: fractal coil geometry + Biot–Savart magnetostatics + moments.
Phase 2: vacuum FDTD (Yee, SSP-RK3) with prescribed coil currents.
Phase 3: self-consistent cold Euler–Maxwell fluid.
"""
module FractalToroid

export Curve, fractal_coil, segments, frames, npoints,
       Box, node, center, cellvol,
       biot_savart, biot_savart_component, biot_savart_yee, vector_potential,
       current_moments, helicity, shell_profile, splat_current!,
       save_png, heatmap_png,
       EX, EY, EZ, BX, BY, BZ, FN, FPX, FPY, FPZ,
       zero_state, em_rhs!, make_sponge, apply_sponge!, field_energy,
       div_B, div_E, ssprk3!, _wrap,
       FluidSim, step!, coupled_rhs!, marder_clean!, kinetic_energy,
       gauss_residual,
       efold_time, powerlaw_slope,
       MRHO, MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI,
       MHDSim, mhd_step!, mhd_rhs!, mhd_kinetic_energy, mhd_magnetic_energy,
       curl_central, grid_moments, azimuthal_spectrum,
       checkpoint_save, checkpoint_load!,
       add_flux_ring!, add_vortex_ring!,
       volume_render

include("geometry.jl")
include("fields.jl")
include("png.jl")
include("yee.jl")
include("fluid.jl")
include("mhd.jl")
include("volren.jl")
include("diagnostics.jl")

end
