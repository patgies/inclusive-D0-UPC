#!/bin/bash
# Roihu modules (the programs need libgsl).

load_roihu_modules() {
	source /usr/share/lmod/lmod/init/bash
	export MODULEPATH=/appl/modulefiles
	module load spack/x86_64/v2026_03/Core/gcc/15.2.0
	module load spack/x86_64/v2026_03/gcc/15.2.0/gsl/2.8
}
