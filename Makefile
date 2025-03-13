#compiler and options
CC = mpiCC
OPTIONS = --std=c++17 -O2 -Wall

#output executable
TARGET = Avcorr.exe

#code files:
SRC = Avcorr.cpp src/CIC.cpp src/Combine.cpp src/Moments.cpp src/Parsfnames.cpp
LIB = lib/Arrvec.cpp lib/Files.cpp lib/Mathfunc.cpp lib/Stat.cpp lib/Strsimsys.cpp

#directories to create if don't exist:
DIRS_add=Data Results Randoms

#to run everything:
all: $(DIRS_add) $(TARGET)

#creating dirs:
$(DIRS_add):
	@mkdir -p $(DIRS_add)

#compilation
$(TARGET): Avcorr.cpp
	$(CC) $(OPTIONS) $(SRC) $(LIB) -o $(TARGET)