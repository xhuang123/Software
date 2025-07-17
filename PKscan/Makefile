pkscan: pkscan.o
	g++ -g -o pkscan pkscan.o

# now define how each object file is a target and list dependencies and how
#  to build that object file if any dependencies change

pkscan.o: pkscan.cpp pkscan.h
	g++ -c pkscan.cpp

clean:
	rm pkscan pkscan.o
