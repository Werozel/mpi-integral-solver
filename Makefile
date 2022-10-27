
TRASH = *.o

TARGET = main.cpp
OUTPUT = main.o

build:
	mpicxx ${TARGET} -o ${OUTPUT}

run:
	mpirun -n ${n_proc} ${OUTPUT} ${eps}

schedule-bluegene:
	mpisubmit.bg -n ${n_proc} ./${OUTPUT}

schedule-polus:
	module load SpectrumMPI/10.1.0
	mpicxx -std=c++11 -O3 -o main.o main.cpp
	mpisubmit.pl -p ${n_proc} main.o ${eps}

clean:
	rm -rf ${TRASH}
