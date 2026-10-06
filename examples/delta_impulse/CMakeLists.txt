# delta_impulse: verification with an INSTANTANEOUS delta impulse.
#
# Reuses the adjoint_forward sources.  Two things differ from the 
# adjoint_forward executable:
#   * src/adjoint_src2_instantaneous.f90 replaces adjoint/adjoint_src2.f90, so
#     the spatial heat kernel is added once at the event time;
#   * INSTANTANEOUS_VERIFICATION_ONLY stops the driver after the unit
#     representer, before the observation-fit/coefficient solve.

set(af    ${PROJECT_SOURCE_DIR}/examples/adjoint_forward)
set(src   ${PROJECT_SOURCE_DIR}/src)
set(src2d ${PROJECT_SOURCE_DIR}/src/2d)

add_library(delta_impulse_f OBJECT
  ${af}/af_common_overlap.f90
  ${af}/af_common_setaux.f90
  ${af}/af_common_qinit.f90

  ${af}/adjoint/adjoint_psi.f90
  ${af}/adjoint/adjoint_fdisc.f90 
  ${af}/adjoint/adjoint_initial.f90
  ${af}/adjoint/adjoint_module.f90
  ${af}/adjoint/adjoint_source.f90
  ${af}/adjoint/adjoint_setprob.f90
  src/adjoint_src2_instantaneous.f90
  ${af}/adjoint/rp/clawpack46_rpn2adv.f
  ${af}/adjoint/rp/clawpack46_rpt2adv.f

  ${af}/forward/forward_psi.f90
  ${af}/forward/forward_initial.f90
  ${af}/forward/forward_module.f90
  ${af}/forward/forward_source.f90
  ${af}/forward/forward_setprob.f90
  ${af}/forward/forward_src2.f90
  ${af}/forward/q_optimal_est.f90
  ${af}/forward/rp/rpn2_forward.f90
  ${af}/forward/rp/rpt2_forward.f90
  ${af}/forward/rp/rpn2_conservative_update.f90
)

# The modules have the same names as adjoint_forward's; keep them apart.
set_target_properties(delta_impulse_f PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
)

# The C/C++ side is identical to adjoint_forward, so its compiled objects
# (adjoint_forward_common, defined in examples/adjoint_forward) are reused.
# Only the driver is compiled here, with the verification-only flag.
add_executable(delta_impulse
  ${af}/adjoint_forward.cpp
  $<TARGET_OBJECTS:adjoint_forward_common>
  $<TARGET_OBJECTS:delta_impulse_f>
)

target_compile_definitions(delta_impulse PRIVATE
  INSTANTANEOUS_VERIFICATION_ONLY=1
)

target_include_directories(delta_impulse PRIVATE
  ${af}
  ${src}
  ${src2d}
)

target_link_libraries(delta_impulse PRIVATE
  smoke2d-4DVar
  FORESTCLAW::CLAWPACK4.6
)
