add_library(
	amplitude
	STATIC
	amplitudelib.cpp
	datafile.cpp
	tools.cpp
	interpolation.cpp
)

add_executable(
    dipole
    inclusive_spectrum.cpp
    main.cpp
	gamma_aa.cpp
    fourier.c
	fragmentation.cpp
	hymnd_grid.cpp
	kk_grid.cpp
	bcfy_grid.cpp
)

target_link_libraries(
	dipole
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
	-static-libstdc++
	-static-libgcc
)


add_executable(
	scan_KK_FF
	scan_KK_FF.cpp
	kk_grid.cpp
	hymnd_grid.cpp
)

target_link_libraries(
	scan_KK_FF
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
	-static-libstdc++
	-static-libgcc
)

set_target_properties(scan_KK_FF PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_executable(
	scan_BCFY_FF
	scan_BCFY_FF.cpp
	bcfy_grid.cpp
	hymnd_grid.cpp
	fragmentation.cpp
)

target_link_libraries(
	scan_BCFY_FF
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
	-static-libstdc++
	-static-libgcc
)

set_target_properties(scan_BCFY_FF PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
