SOURCE	= fft99_modified.f90 parameters.f90 compute_divergence_tendency.f90 compute_vorticity_tendency.f90  init.f90 \
legt_i.f90 save_output.f90 compute_geopotential_tendency.f90 convert_vordiv2uv.f90 convert_uv2vordiv.f90 fft_d.f90 \
j_index2.f90 main_sswm.f90 compute_kinetic_energy.f90 compute_ke_spectrum.f90 d_legpol.f90 fft_i.f90 legt_d.f90 

FC = gfortran
#FFLAGS = -fdefault-real-8	
#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o sswe.exe ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) *~

fort:
	 \rm -f  fort.* *~
