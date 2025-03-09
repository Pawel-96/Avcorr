#compiler and options
CC = mpiCC
OPTIONS = --std=c++17 -O2 -Wall

#output executable
TARGET = Avcorr.exe

#code files:
SRC = Avcorr.cpp src/CIC.cpp src/Combine.cpp src/Moments.cpp src/Parsfnames.cpp
LIB = lib/Arrvec.cpp lib/Files.cpp lib/Mathfunc.cpp lib/Stat.cpp lib/Strsimsys.cpp

#compilation
$(TARGET): Avcorr.cpp
	$(CC) $(OPTIONS) $(SRC) $(LIB) -o $(TARGET)