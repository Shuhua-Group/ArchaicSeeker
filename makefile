all: ArchaicSeeker

ArchaicSeeker: ArchData.o ArchSeeker.o ArchPar.o
	g++ ArchData.o ArchSeeker.o ArchPar.o -o ArchaicSeeker -lz -lm

ArchSeeker.o: ArchSeeker.cpp
	g++ -c ArchSeeker.cpp

ArchData.o: ArchData.cpp ArchData.h ArchTools.h
	g++ -c ArchData.cpp -lz

ArchPar.o: ArchPar.cpp ArchPar.h ArchTools.h
	g++ -c ArchPar.cpp

clean:
	rm -f *.o ArchaicSeeker

.PHONY: all clean
