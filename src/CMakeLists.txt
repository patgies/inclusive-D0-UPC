
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
	lhapdf_grid.cpp
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
