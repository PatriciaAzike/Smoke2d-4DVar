# --------------------------------------------
# List of available applications and examples
#
# Note : Any files below which are uncommented
# need to be available at compile time, even 
# the configuration flag is not set.  
# --------------------------------------------


cmake_minimum_required(VERSION 3.13...3.20)

project(Smoke2d-4DVar LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

# # -- check if this being included as a subdirectory
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ForestClawExternalProject)

# --- get installed ForestClaw package

# Set to True if you want to download ForestClaw and build from scratch
if(NOT external_build)
    # 'clawpack' here refers to the "-Dclawpack=on" cmake flag
    # This uses -DFORESTCLAW_ROOT to find ForestClaw installation
    find_package(FORESTCLAW REQUIRED COMPONENTS clawpack)
else()
    # --- fall back to external project build.  ForestClaw will be cloned in the
    # build directory.  Requires internet access to GitHub. 
    include(FetchContent)
    FetchContent_Declare(forestclaw
                    GIT_REPOSITORY https://www.github.com/ForestClaw/forestclaw.git
                    GIT_TAG develop)
    set(mpi TRUE)
    set(FCLAW_ENABLE_MPI TRUE)
    FetchContent_MakeAvailable(forestclaw)
endif()
 
add_subdirectory(src)
add_subdirectory(examples)
