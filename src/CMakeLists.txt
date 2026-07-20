
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
	# Statically link the C++ runtime so the binary doesn't depend on
	# whatever libstdc++.so.6 happens to be installed on a given compute
	# node -- login/build nodes often have a newer GCC module loaded than
	# compute nodes get by default, causing GLIBCXX_* "version not found"
	# errors at runtime otherwise.
	-static-libstdc++
	-static-libgcc
)
