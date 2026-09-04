
# QCDNum (used by kk_grid.cpp to DGLAP-evolve the Kniehl & Kramer
# fragmentation function). Vendored under QCDnum/ (header + static libs
# only) so this project builds standalone, without needing a QCDNum
# source/build tree present on the machine. Override with -DQCDNUM_DIR=...
# to point at a different build instead.
set(QCDNUM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../QCDnum" CACHE PATH
	"Path to a QCDNum install (containing include/QCDNUM/QCDNUM.h, lib/*.a)")

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
)

target_include_directories(
	dipole
	PRIVATE
	${QCDNUM_DIR}/include
)

target_link_libraries(
	dipole
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
	${QCDNUM_DIR}/lib/libQCDNUMCCWrap.a
	${QCDNUM_DIR}/lib/libQCDNUM.a
	gfortran
	-static-libstdc++
	-static-libgcc
)
