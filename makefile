OBJ= precision.o  constants.o systems.o channels.o mesh.o pot.o input.o   coul90.o gauss.o interpolation.o  lagrange_mesh.o scattwf.o scattering.o 


# Laptop
 LIB =-L  /opt/local/lib -llapack
#LIB = -L ../lapack-3.5.0 -lrefblas -llapack  
FC = gfortran 
F90 = gfortran
FFLAGS = -O2 -Wtabs   -ffixed-line-length-0  

.SUFFIXES: .F90 .f90 .f95

all: scattering


scattering:  $(OBJ)
	$(FC) -o scattering $(FFLAGS) $(OBJ) $(LIB)




.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f95.o          :
	$(F90) $(FFLAGS) -c $<

.F.o          :
	$(F90) $(FFLAGS) -c $<

.f.o          :
	$(F90) $(FFLAGS) -c $<


clean: 
	rm -f upot $(objectsbound) $(objectsscatt) *.mod core *.o scattering


