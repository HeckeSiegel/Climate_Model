#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f90 )

SET_SOURCE_FILES_PROPERTIES(rrtmg_sw_k_g.f90 PROPERTIES COMPILE_FLAGS -O0)

add_library( fpda_rrtm_sw STATIC ${files} )
TARGET_LINK_LIBRARIES( fpda_rrtm_sw )

add_executable( fpda_example_sw fpda_example.F90 )
target_link_libraries( fpda_example_sw fpda_rrtm_sw )

install(TARGETS fpda_rrtm_sw   DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
