PROGRAM = burster
LOADLIBES = -lm
CXX.SRCS = ran3.c main.c
CC = gcc
LDFLAGS = -O3

CXX = -O3
CXXFLAGS = -O3 -Wall $(DEFINES)

OBJS = $(CXX.SRCS:.cc=.o) 

$(PROGRAM) : $(OBJS)
	$(CC) $(LDFLAGS) $(DEFINES) $(OBJS) $(LOADLIBES) -o $(PROGRAM).x

clean:
	/bin/rm -f *.o $(PROGRAM).x *~

depend:
	$(CXX) -MM $(CXX.SRCS)

### DEPENDENCIES ON .h FILES GO HERE

$(OBJS) : main.h
