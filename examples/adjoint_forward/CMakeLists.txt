set(src ${PROJECT_SOURCE_DIR}/src)
set(src2d ${PROJECT_SOURCE_DIR}/src/2d)

add_library(adjoint_forward_f OBJECT 
  af_common_overlap.f90
  af_common_setaux.f90
  af_common_qinit.f90

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
  user_run.c

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
 smoke2d-4DVar
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
 smoke2d-4DVar
 FORESTCLAW::CLAWPACK4.6
)

# The instantaneous delta-impulse verification is its own example:
# examples/delta_impulse (it reuses the sources above).
