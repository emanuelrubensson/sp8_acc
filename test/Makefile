TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/source)
LIBS=-llapack
all: test_homotopy_solver test_sp8_evaluation_scalar

test_homotopy_solver: test_homotopy_solver.cc
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_sp8_evaluation_scalar: test_sp8_evaluation_scalar.cc 
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_eval_sp8: test_eval_sp8.cc
	$(CXX) $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

check: test_homotopy_solver test_sp8_evaluation_scalar
	./test_homotopy_solver > test_homotopy_solver_output.txt && grep -q 'finished OK!' test_homotopy_solver_output.txt || exit 1;
	rm test_homotopy_solver_output.txt
	printf "test_homotopy_solver finished OK!\n"
	./test_sp8_evaluation_scalar > test_sp8_evaluation_scalar_output.txt && grep -q 'finished OK!' test_sp8_evaluation_scalar_output.txt || exit 1;
	rm test_sp8_evaluation_scalar_output.txt
	printf "test_sp8_evaluation_scalar finished OK!\n"
clean:
	rm test_homotopy_solver test_sp8_evaluation_scalar test_eval_sp8
