TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/sp8cc)
LIBS=-llapack -lblas
all: test_homotopy_solver test_sp8_evaluation_scalar test_poly8_evaluation_matrix test_poly8_evaluation_matrix_symmetric test_get_sp8_params_max_gap

test_homotopy_solver: test_homotopy_solver.cc
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_sp8_evaluation_scalar: test_sp8_evaluation_scalar.cc 
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_poly8_evaluation_matrix: test_poly8_evaluation_matrix.cc $(TOPDIR)/sp8cc/poly8_eval_low_mem.h
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_poly8_evaluation_matrix_symmetric: test_poly8_evaluation_matrix_symmetric.cc
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)
test_get_sp8_params_max_gap: test_get_sp8_params_max_gap.cc
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_eval_sp8: test_eval_sp8.cc
	$(CXX) $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

check: all
	./test_homotopy_solver > test_homotopy_solver_output.txt && grep -q 'finished OK!' test_homotopy_solver_output.txt || exit 1;
	rm test_homotopy_solver_output.txt
	printf "test_homotopy_solver finished OK!\n"
	./test_sp8_evaluation_scalar > test_sp8_evaluation_scalar_output.txt && grep -q 'finished OK!' test_sp8_evaluation_scalar_output.txt || exit 1;
	rm test_sp8_evaluation_scalar_output.txt
	printf "test_sp8_evaluation_scalar finished OK!\n"
	./test_poly8_evaluation_matrix > test_poly8_evaluation_matrix_output.txt && grep -q 'finished OK!' test_poly8_evaluation_matrix_output.txt || exit 1;
	rm test_poly8_evaluation_matrix_output.txt
	printf "test_poly8_evaluation_matrix finished OK!\n"
	./test_poly8_evaluation_matrix_symmetric > test_poly8_evaluation_matrix_symmetric_output.txt && grep -q 'finished OK!' test_poly8_evaluation_matrix_symmetric_output.txt || exit 1;
	rm test_poly8_evaluation_matrix_symmetric_output.txt
	printf "test_poly8_evaluation_matrix_symmetric finished OK!\n"
	./test_get_sp8_params_max_gap > test_get_sp8_params_max_gap_output.txt && grep -q 'finished OK!' test_get_sp8_params_max_gap_output.txt || exit 1;
	rm test_get_sp8_params_max_gap_output.txt
	printf "test_get_sp8_params_max_gap finished OK!\n"
clean:
	rm test_homotopy_solver test_sp8_evaluation_scalar test_poly8_evaluation_matrix test_poly8_evaluation_matrix_symmetric test_eval_sp8 test_get_sp8_params_max_gap
