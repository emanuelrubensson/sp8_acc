TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/source)
LIBS=-llapack
test_solver_43: test_solver_43.cc
	$(CXX) $(SP8_ACC_INCL) -o $@ $^ $(LIBS) 
