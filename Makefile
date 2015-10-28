CXX       = g++
CXXFLAGS  = -O2 -std=c++11

pixel_tracing: pixel_tracing.cpp
	$(CXX) -ggdb3 -o $@ pixel_tracing.cpp -lSDL

.depend:
	$(CXX) $(CXXFLAGS) $(INC) -MM *.cpp > .depend

-include .depend
