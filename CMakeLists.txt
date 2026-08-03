cmake_minimum_required(VERSION 3.16)

project(DG LANGUAGES Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_executable(IFS_01_preprocessing
    src/01_preprocess/01_module.f90
    src/01_preprocess/01_read.f90
    src/01_preprocess/01_mesh_dg.f90
    src/01_preprocess/01_post.f90
    src/01_preprocess/01_deallocate.f90
    src/01_preprocess/01_main.f90
)
set_target_properties(IFS_01_preprocessing PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/01_preprocess
)
target_include_directories(IFS_01_preprocessing PRIVATE ${CMAKE_BINARY_DIR}/mod/01_preprocess)

add_executable(IFS_02_solver_dg
    src/02_solver_dg/02_module.f90
    src/02_solver_dg/02_read.f90
    src/02_solver_dg/02_basis.f90
    src/02_solver_dg/02_init.f90
    src/02_solver_dg/02_conv.f90
    src/02_solver_dg/02_diff.f90
    src/02_solver_dg/02_time.f90
    src/02_solver_dg/02_deallocate.f90
    src/02_solver_dg/02_main.f90
)
set_target_properties(IFS_02_solver_dg PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/02_solver_dg
)
target_include_directories(IFS_02_solver_dg PRIVATE ${CMAKE_BINARY_DIR}/mod/02_solver_dg)
