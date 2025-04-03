#compiler and options
CC = mpiCC
OPTIONS = --std=c++17 -lhdf5_cpp -O2 -Wall
H5flag=`pkg-config --cflags hdf5-serial`
H5link=`pkg-config --libs hdf5-serial`

#output executable
TARGET = Avcorr.exe

#code files:
SRC = Avcorr.cpp src/CIC.cpp src/Combine.cpp src/Moments.cpp src/Parsfnames.cpp
LIB = lib/Arrvec.cpp lib/Files.cpp lib/Mathfunc.cpp lib/Stat.cpp lib/Strsimsys.cpp
H5FILE = lib/Files_HDF5.cpp

#directories to create if don't exist:
DIRS_add=Data Results Randoms

#to run everything:
all: $(DIRS_add) $(TARGET)

#creating dirs:
$(DIRS_add):
	@mkdir -p $(DIRS_add)

#compilation
$(TARGET): Avcorr.cpp
	$(CC) $(H5flag) $(H5FILE) $(H5link) $(OPTIONS) $(SRC) $(LIB) -o $(TARGET)
