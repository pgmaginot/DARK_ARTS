CXXFLAGS = -Wall -std=c++11
LIBS = 
SRCS = tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp Input_Reader.cc input_test_main.cc
INCLUDES =  tinystr.h tinyxml.h Input_Reader.h Inputs_Allowed.h
OBJS = ${SRCS:.cpp=.o}
TARG = inpTest
CXX = g++

all:
	@$(MAKE) CXXFLAGS="$(CXXFLAGS) -O3" TARG="$(TARG)" $(TARG)

debug:
	@$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" TARG="$(TARG)" $(TARG)

debug_opt:
	@$(MAKE) CXXFLAGS="$(CXXFLAGS) -g -O3" TARG="$(TARG)" $(TARG)

clean:
	/bin/rm -rf $(OBJS) $(TARG) core
	
$(TARG): $(INCLUDES) $(OBJS)
	$(CXX) -o $(TARG) $(CXXFLAGS) $(OBJS) $(LIBS)