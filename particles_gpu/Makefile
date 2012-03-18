# Load CUDA using the following command
# module load cuda

#
# Dirac
#
CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch sm_20
NVCCFLAGS = -O3 -arch sm_20
LIBS = -lm

TARGETS = serial gpu

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
gpu: gpu.o common.o
	$(CC) -o $@ $(NVCCLIBS) gpu.o common.o

serial.o: serial.cu common.h
	$(CC) -c $(CFLAGS) serial.cu
gpu.o: gpu.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) common.cu

clean:
	rm -f *.o $(TARGETS)
