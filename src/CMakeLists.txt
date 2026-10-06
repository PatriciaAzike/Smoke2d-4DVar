
add_library(smoke2d_f OBJECT 
  smoke2d_module.f90
  smoke2d_gauges_fort.f90
  2d/transport2d_psi.f90
  2d/transport2d_velocity.f90
  2d/transport2d_setaux_velocity.f90
  #2d/rpn2_forward.f90
  #2d/rpt2_forward.f90
  #2d/rpn2_conservative_update.f90
  #2d/rpn2_adjoint.f90
  #2d/rpt2_adjoint.f90
  2d/src2_adjoint.f90

  fortran_source/clawpatch46_average.f90
  fortran_source/clawpatch46_copy.f90
  fortran_source/clawpatch46_fort_tag4coarsening.f90
  fortran_source/clawpatch46_fort_tag4refinement.f90
  fortran_source/clawpatch46_interpolate.f
  fortran_source/clawpatch46_local_ghost_pack.f
  fortran_source/clawpatch46_timeinterp.f
)

# Needed for modules
target_include_directories(smoke2d_f PUBLIC 
     ${CMAKE_CURRENT_BINARY_DIR}
     ${CMAKE_CURRENT_BINARY_DIR/2d}
)

add_library(smoke2d-4DVar
  smoke2d_options.c
  smoke2d_gauges.c
  2d/transport2d_setaux.c
  2d/transport2d_b4step2.c
  $<TARGET_OBJECTS:smoke2d_f>  
)

target_link_libraries(smoke2d-4DVar PRIVATE 
 FORESTCLAW::CLAWPACK4.6 
)


