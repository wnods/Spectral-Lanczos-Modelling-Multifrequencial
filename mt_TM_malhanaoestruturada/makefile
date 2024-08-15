

#link = -L /home/zoraida/pacotes_pardiso -lpardiso600-GNU800-X86-64 \
#-L /home/zoraida/pacotes_pardiso/lapack-3.9.0/BLAS \
#-l /home/zoraida/pacotes_pardiso/lapack-3.9.0/SRC -fopenmp -lpthread -lm

#link = -L/home/zoraida/pacotes_pardiso/ -lpardiso600-GNU800-X86-64 \
#-L/home/zoraida/pacotes_pardiso/lapack-3.9.0/ -llapack -lrefblas -fopenmp -lpthread -lm

link = -lmkl_rt

run:
	gfortran modulo_varglob.f90 mod_arrays_pds_cmplx.f90 modulo_femmttm.f90 principal.f90 -o mt2D.x $(link)
expt:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/zoraida/pacotes_pardiso/
	export PARDISOLICMESSAGE=1
clean:
	rm *.mod mt2D.x
