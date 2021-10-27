TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/source)
LIBS=-llapack
test_homotopy_solver: test_homotopy_solver.cc
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_eval_sp8: test_eval_sp8.cc
	$(CXX) $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

check: test_homotopy_solver
	./test_homotopy_solver > test_homotopy_solver_output.txt && grep -q 'finished OK!' test_homotopy_solver_output.txt || exit 1;
	rm test_homotopy_solver_output.txt
	printf "test_homotopy_solver finished OK!\n"
clean:
	rm test_homotopy_solver test_eval_sp8
