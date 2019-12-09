#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

#add_executable( rrtm_example_lw_FJ fpda_example_lw_FJ.c )
#target_link_libraries( rrtm_example_lw_FJ m fpda_rrtm_lw gnuplot_i )

add_executable( rrtm_example_lw_BM fpda_example_lw_BM.c )
target_link_libraries( rrtm_example_lw_BM m fpda_rrtm_lw gnuplot_i )

#add_executable( rrtm_example_sw_FJ fpda_example_sw_FJ.c )
#target_link_libraries( rrtm_example_sw_FJ m fpda_rrtm_sw gnuplot_i )

add_executable( rrtm_example_sw_BM fpda_example_sw_BM.c )
target_link_libraries( rrtm_example_sw_BM m fpda_rrtm_sw gnuplot_i )


#install(TARGETS rrtm_example_lw_FJ DESTINATION bin)
install(TARGETS rrtm_example_lw_BM DESTINATION bin)
#install(TARGETS rrtm_example_sw_FJ DESTINATION bin)
install(TARGETS rrtm_example_sw_BM DESTINATION bin)
