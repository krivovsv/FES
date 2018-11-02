npvrcop: npvrcop.f npvrcoplib.f
	gfortran -O3 -ffast-math -march=native -funroll-loops  -o npvrcop npvrcop.f npvrcoplib.f -llapack -fopenmp
