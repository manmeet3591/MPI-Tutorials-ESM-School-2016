CC = mpicc
FC = mpiifort
all:
	$(FC) dotproduct.f90 -o dotproduct.x

extra:
	$(CC) extra_dotproduct.c -o extra_dotproduct.x

serial:	
#	gcc serial_dotproduct.c -o serial_dotproduct.x
	ifort serial_dotproduct.f90 -o serial_dotproduct.x
