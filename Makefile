CXX       = g++

G2O_DIR=/home/jd/Downloads/g2o-master

all: pixel_tracing test_g2o improc match_g2o viz_cam match_solver plot_epipole check_epipole linear_triangulate abs_traj_err

pixel_tracing: pixel_tracing.cpp pose_est.cpp geometry.h
	$(CXX) -ggdb3 -o $@ pixel_tracing.cpp pose_est.cpp -lSDL

test_g2o: test_g2o.cpp
	$(CXX) -std=c++11 -ggdb3 -I$(G2O_DIR) -I/usr/local/include/eigen3 -L$(G2O_DIR)/lib -Wl,-R$(G2O_DIR)/lib -o $@ $< -lg2o_core -lg2o_solver_dense -lg2o_stuff

improc: improc.c
	gcc -Wall -std=c99 -ggdb3 -o $@ $< -llapack

match_g2o: match_g2o.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -I$(G2O_DIR) -I/usr/local/include/eigen3 -L$(G2O_DIR)/lib -Wl,-R$(G2O_DIR)/lib -o $@ $< -lg2o_core -lg2o_solver_dense -lg2o_stuff

match_solver: match_solver.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -I/usr/local/include/eigen3 -o $@ $<

linear_triangulate: linear_triangulate.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -o $@ $< -lopencv_core

abs_traj_err: abs_traj_err.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -o $@ $< -I/usr/local/include/eigen3

viz_cam: viz_cam.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -o $@ $< -lSDL -lGL -lGLU

plot_epipole: plot_epipole.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -o $@ $< -lSDL -lGL -lGLU

check_epipole: check_epipole.cpp
	$(CXX) -Wall -std=c++11 -ggdb3 -o $@ $<
