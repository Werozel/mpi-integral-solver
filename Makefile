
TRASH = *.o

TARGET = main.cpp
OUTPUT = main.o

build:
	mpicxx ${TARGET} -o ${OUTPUT}

run:
	mpirun -n ${n_proc} ${OUTPUT}

schedule:
	mpisubmit.bg -n ${n_proc} ./${OUTPUT}

clean:
	rm -rf ${TRASH}
