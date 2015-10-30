CXX       = g++
CXXFLAGS  = -O2 -std=c++11

pixel_tracing: pixel_tracing.cpp pose_est.cpp geometry.h
	$(CXX) -ggdb3 -o $@ pixel_tracing.cpp pose_est.cpp -lSDL

.depend:
	$(CXX) $(CXXFLAGS) $(INC) -MM *.cpp > .depend

-include .depend
