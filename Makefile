all: npvrcop npvevop

npvrcop: npvrcop.f npvrcoplib.f
	gfortran -O3 -ffast-math -march=native -funroll-loops  -o npvrcop npvrcop.f npvrcoplib.f -llapack -fopenmp
npvevop: npvevop.f npvevoplib.f
	gfortran -O3 -ffast-math -march=native -funroll-loops  -o npvevop npvevop.f npvevoplib.f -llapack -fopenmp
#	 -fbounds-check -g
