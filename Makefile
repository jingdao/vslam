CXX       = g++

G2O_DIR=/home/jd/Downloads/g2o-master

all: match_g2o

pixel_tracing: pixel_tracing.cpp pose_est.cpp geometry.h
	$(CXX) -ggdb3 -o $@ pixel_tracing.cpp pose_est.cpp -lSDL

test_g2o: test_g2o.cpp
	$(CXX) -std=c++11 -ggdb3 -I$(G2O_DIR) -I/usr/local/include/eigen3 -L$(G2O_DIR)/lib -Wl,-R$(G2O_DIR)/lib -o $@ $< -lg2o_core -lg2o_solver_dense -lg2o_stuff

improc: improc.c
	gcc -Wall -std=c99 -ggdb3 -o $@ $< -llapack

match_g2o: match_g2o.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -I$(G2O_DIR) -I/usr/local/include/eigen3 -L$(G2O_DIR)/lib -Wl,-R$(G2O_DIR)/lib -o $@ $< -lg2o_core -lg2o_solver_dense -lg2o_stuff

