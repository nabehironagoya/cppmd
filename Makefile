CC = gcc
CXX = g++
CFLAGS =  
CXXFLAGS =  -g -O3 -std=gnu++11 -Wall #-DDEBUG  # -I /opt/apl/gromacs2019.2/include/ -I /misc/home/tnagai/src/cxx-prettyprint
LOADLIBES =  
#LDLIBS_BOOST =  -lboost_program_options
#LDLIBS_GROMACS =  -L /opt/apl/gromacs2019.2/lib64/ -lgromacs 

OBJS1 =   read_gro.o 
OBJS2 =   string.o 
OBJS3 =   math.o 
OBJS4 =   mdrun.o 
OBJ_NAMES = $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)

EXE_NAME1 = cppmd
EXE_NAMES = $(EXE_NAME1)  

all: $(EXE_NAMES)

$(EXE_NAME1): $(OBJ_NAMES)
	$(CXX) -o $(EXE_NAME1)  $(OBJ_NAMES)  # $(LOADLIBES)  $(LDLIBS_BOOST)


.PHONY: clean uninstall rebuild 
clean:
	    rm -f $(OBJ_NAMES)

uninstall:
	    rm -f $(OBJ_NAMES) $(EXE_NAMES)

rebuild: uninstall $(EXE_NAMES) 


