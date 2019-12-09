#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_library( gnuplot_i STATIC gnuplot_i.c )
TARGET_LINK_LIBRARIES( gnuplot_i )

install(TARGETS gnuplot_i DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")

