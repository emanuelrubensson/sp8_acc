TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/source)
LIBS=-llapack
test_solver_43: test_solver_43.cc
	$(CXX) $(SP8_ACC_INCL) -o $@ $^ $(LIBS) 

check: test_solver_43
	./test_solver_43 > test_solver_43_output.txt && grep -q 'finished OK!' test_solver_43_output.txt || exit 1;
	rm test_solver_43_output.txt
	printf "test_solver_43 finished OK!\n"
clean:
	rm test_solver_43
