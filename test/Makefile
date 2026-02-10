TOPDIR=..
SP8_ACC_INCL=-I$(realpath $(TOPDIR)/sp8cc)
SP8_ACC_HEADER_FILES := $(wildcard $(TOPDIR)/sp8cc/*.h)
LIBS=-llapack -lblas
all: test_homotopy_solver test_sp8_evaluation_scalar test_poly8_evaluation_matrix test_poly8_evaluation_matrix_symmetric test_get_sp8_params_max_gap

test_homotopy_solver: test_homotopy_solver.cc \
	$(SP8_ACC_HEADER_FILES)
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_sp8_evaluation_scalar: test_sp8_evaluation_scalar.cc \
	$(SP8_ACC_HEADER_FILES)
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_poly8_evaluation_matrix: test_poly8_evaluation_matrix.cc \
	$(SP8_ACC_HEADER_FILES)
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_poly8_evaluation_matrix_symmetric: test_poly8_evaluation_matrix_symmetric.cc \
	$(SP8_ACC_HEADER_FILES)
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_get_sp8_params_max_gap: test_get_sp8_params_max_gap.cc \
	$(SP8_ACC_HEADER_FILES)
	$(CXX) -std=c++17 $(SP8_ACC_INCL) -o $@ $^ $(LIBS)

test_eval_sp8: test_eval_sp8.cc \
	$(SP8_ACC_HEADER_FILES)
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
	rm -f test_homotopy_solver \
	test_sp8_evaluation_scalar \
	test_poly8_evaluation_matrix \
	test_poly8_evaluation_matrix_symmetric \
	test_eval_sp8 \
	test_get_sp8_params_max_gap