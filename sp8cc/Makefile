
HEADER_FILES=Complex_step_solver.h sp8_acc.h
LIBS=-llapack
sp8cc.so: $(HEADER_FILES) cwrappers.cc
	$(CXX) -std=c++17 -fPIC -shared -o $@ cwrappers.cc $(LIBS)
