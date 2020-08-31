CC = g++
MyO = -I. -I./ranlxd -I$(CPATH) -DDSFMT_MEXP=19937 -DUSE_RAND55 -D__PARALLEL_OMP__ -g -O3 -std=c++0x -fopenmp
MyL = -L. -L$(LIBRARY_PATH) -lm -lfftw3 -lfftw3_threads -fopenmp

Quenched: libprlgt.a nspt.o Quenched.o *.h 
	  $(CC) -o Quenched nspt.o Quenched.o libprlgt.a $(MyL)

libprlgt.a: MyMath.o MyRand.o MyQCD.o QCDpt.o ranlxd/ranlxd.o *.h
	  ar cru libprlgt.a MyMath.o MyRand.o MyQCD.o QCDpt.o ranlxd/ranlxd.o

nspt.o:
	  $(CC) $(MyO) -c -o nspt.o nspt.cc

Quenched.o: 
	  $(CC) $(MyO) -c -o Quenched.o Quenched.cc

MyMath.o: 
	  $(CC) $(MyO) -c -o MyMath.o MyMath.cc

MyRand.o: 
	  $(CC) $(MyO) -c -o MyRand.o MyRand.cc

MyQCD.o: 
	  $(CC) $(MyO) -c -o MyQCD.o MyQCD.cc

QCDpt.o: 
	  $(CC) $(MyO) -c -o QCDpt.o QCDpt.cc

ranlxd/ranlxd.o:
	  $(CC) $(MyO) -c -o ranlxd/ranlxd.o ranlxd/ranlxd.cc


.PHONY : clean

clean:
	  rm -f Quenched *.o *.a */*.o

