SRCS = main.cpp ply_io.c icp.cpp
OUTPUT = ICP_stitch
LIBS = -lboost_system -lboost_filesystem -lopencv_core

all:
	g++ -std=c++11 $(SRCS) -o $(OUTPUT) -Wno-write-strings $(LIBS)

clean: 
	rm -f $(OUTPUT) 
