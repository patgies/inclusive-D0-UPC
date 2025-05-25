
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
        function.cpp
        integral_parton.cpp
        integral_D0.cpp
        main.cpp
        fourier.c
        #libFcompat.c
)	

target_link_libraries(
	dipole
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
)
