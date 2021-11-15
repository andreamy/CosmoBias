cfastcov_dir := ../cosmolike_core/2dfftlog/
cfastcov := $(cfastcov_dir)twobessel.c $(cfastcov_dir)utils.c $(cfastcov_dir)utils_complex.c
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99


cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

################### CLASS SECTION, aka Julien Lesgourges CLASS##################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Julien's class must lie in a folder "class" WITHOUT version number
# curently it is version 3.0

#Directory in which Julien Lesgourges et al's CLASS lies
J:= ../cosmolike_core/class/
B:= $(J)build/

INC=-I$(J)include -I$(J)cpp

#FOR CLASS v3.0
#because it allows the CLASS developers to keep their relative path names.
OBJ=$(B)arrays.o $(B)background.o $(B)common.o $(B)dei_rkck.o \
$(B)distortions.o $(B)energy_injection.o $(B)evolver_ndf15.o \
$(B)evolver_rkck.o $(B)fourier.o $(B)growTable.o $(B)harmonic.o \
$(B)helium.o $(B)history.o  $(B)hydrogen.o $(B)hyperspherical.o \
$(B)hyrectools.o $(B)injection.o $(B)input.o $(B)lensing.o \
$(B)noninjection.o $(B)output.o $(B)parser.o $(B)perturbations.o  \
$(B)primordial.o $(B)quadrature.o $(B)sparse.o $(B)thermodynamics.o \
$(B)transfer.o $(B)trigonometric_integrals.o $(B)wrap_hyrec.o $(B)wrap_recfast.o

#If you have made changes to Julien's class, recompile it:
#This compiles the stand-alone version of class
#class0:
#	cd $(J) && make class

#This compiles the class engine, such that external codes can communicate
#with class
#class1:
#	cd $(J)cpp && gcc -O2 -fopenmp $(INC) -c Engine.cc -o Engine.o \
#	&& gcc -O2 -fopenmp $(INC) -c ClassEngine.cc -o ClassEngine.o

#clean:
#	rm *.o && rm $(E)*.o && rm $(S)*.o

covs:
	make cov_fft
	make cov_curved

cov_fft:
#	gcc -o cov_flat_fft compute_covariances_real_flat_fft.c $(cfastcov) $(opt_home)
	gcc -c compute_covariances_real_flat_fft.c -fopenmp $(INC) $(cfastcov) $(opt_home)
	gcc -o cov_flat_fft compute_covariances_real_flat_fft.o $(OBJ) -fopenmp $(cfastcov) $(opt_home)

cov_curved:
	gcc -o cov compute_covariances_real.c $(cfftlog) $(opt_home)

response:
	gcc -o cov_response compute_covariances_real_WL_response.c  $(opt_home) -lpython2.7


opt_nersc := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I ${GSL_DIR}/include  -L ${GSL_DIR}/lib -I ${FFTW_INC} -L ${FFTW_DIR} -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99

nersc_covs:
	make nersc_cov_fft
	make nersc_cov_curved

nersc_cov_fft:
	gcc -o cov_flat_fft compute_covariances_real_flat_fft.c $(cfastcov) $(opt_nersc)

nersc_cov_curved:
	gcc -o cov compute_covariances_real.c $(cfftlog) $(opt_nersc)
