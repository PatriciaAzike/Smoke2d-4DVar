set(src ${PROJECT_SOURCE_DIR}/src)
set(src2d ${PROJECT_SOURCE_DIR}/src/2d)

add_library(adjoint_forward_f OBJECT 
  af_common_overlap.f90
  af_common_setaux.f90
  af_common_qinit.f90

  fortran_source/clawpatch46_average.f90
  fortran_source/clawpatch46_copy.f90
  fortran_source/clawpatch46_fort_tag4coarsening.f90
  fortran_source/clawpatch46_fort_tag4refinement.f90
  fortran_source/clawpatch46_interpolate.f
  fortran_source/clawpatch46_local_ghost_pack.f
  fortran_source/clawpatch46_timeinterp.f

  adjoint/adjoint_psi.f90 
  adjoint/adjoint_fdisc.f90 
  adjoint/adjoint_initial.f90
  adjoint/adjoint_module.f90
  adjoint/adjoint_source.f90
  adjoint/adjoint_setprob.f90 
  adjoint/adjoint_src2.f90 
  adjoint/rp/clawpack46_rpn2adv.f 
  adjoint/rp/clawpack46_rpt2adv.f 

  forward/forward_psi.f90 
  forward/forward_initial.f90
  forward/forward_module.f90
  forward/forward_source.f90
  forward/forward_setprob.f90 
  forward/forward_src2.f90 
  forward/q_optimal_est.f90
  forward/rp/rpn2_forward.f90
  forward/rp/rpt2_forward.f90 
  forward/rp/rpn2_conservative_update.f90  
)

add_library(adjoint_forward_common OBJECT 
  af_overlap_patch.c
  user_run_2.c

  adjoint/adjoint.cpp
  adjoint/adjoint_options.c 
  adjoint/adjoint_user.cpp 

  forward/forward.cpp
  forward/forward_user.cpp 
  forward/forward_options.c
  forward/model_options.c
)

target_include_directories(adjoint_forward_common PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${src}
  ${src2d}
)

target_link_libraries(adjoint_forward_common PRIVATE 
 smoke2d-4dVar
 FORESTCLAW::CLAWPACK4.6
)

add_executable(adjoint_forward 
  adjoint_forward.cpp
  $<TARGET_OBJECTS:adjoint_forward_common>
  $<TARGET_OBJECTS:adjoint_forward_f>
)

target_include_directories(adjoint_forward PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${src}
  ${src2d}
)

target_link_libraries(adjoint_forward PRIVATE 
 #smoke3d
 FORESTCLAW::CLAWPACK4.6
)

# Isolated verification executable for a point observation applied as a
# discrete event in time.  It deliberately does not change the production
# adjoint_forward target, whose source is regularized in both space and time.
add_library(adjoint_forward_instantaneous_f OBJECT
  af_common_overlap.f90
  af_common_setaux.f90
  af_common_qinit.f90

  fortran_source/clawpatch46_average.f90
  fortran_source/clawpatch46_copy.f90
  fortran_source/clawpatch46_fort_tag4coarsening.f90
  fortran_source/clawpatch46_fort_tag4refinement.f90
  fortran_source/clawpatch46_interpolate.f
  fortran_source/clawpatch46_local_ghost_pack.f
  fortran_source/clawpatch46_timeinterp.f

  adjoint/adjoint_psi.f90
  adjoint/adjoint_fdisc.f90
  adjoint/adjoint_initial.f90
  adjoint/adjoint_module.f90
  adjoint/adjoint_source.f90
  adjoint/adjoint_setprob.f90
  delta_impulse_instantaneous/src/adjoint_src2_instantaneous.f90
  adjoint/rp/clawpack46_rpn2adv.f
  adjoint/rp/clawpack46_rpt2adv.f

  forward/forward_psi.f90
  forward/forward_initial.f90
  forward/forward_module.f90
  forward/forward_source.f90
  forward/forward_setprob.f90
  forward/forward_src2.f90
  forward/q_optimal_est.f90
  forward/rp/rpn2_forward.f90
  forward/rp/rpt2_forward.f90
  forward/rp/rpn2_conservative_update.f90
)

# Both Fortran object libraries compile modules with the same names.  Keep the
# alternate module files out of the production target's module directory 
set_target_properties(adjoint_forward_instantaneous_f PROPERTIES
  Fortran_MODULE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/instantaneous_fortran_modules
)

add_executable(adjoint_forward_instantaneous
  adjoint_forward.cpp
  $<TARGET_OBJECTS:adjoint_forward_common>
  $<TARGET_OBJECTS:adjoint_forward_instantaneous_f>
)

# This target verifies one unit representer only.  Do not continue into the
# observation-fit/coefficient solve after the adjoint-forward propagation.
target_compile_definitions(adjoint_forward_instantaneous PRIVATE
  INSTANTANEOUS_VERIFICATION_ONLY=1
)

target_include_directories(adjoint_forward_instantaneous PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${src}
  ${src2d}
)

target_link_libraries(adjoint_forward_instantaneous PRIVATE
  FORESTCLAW::CLAWPACK4.6
)
