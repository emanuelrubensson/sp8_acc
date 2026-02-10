
HEADER_FILES=Complex_step_solver.h sp8_acc.h Matrix_proxy.h sastre_poly8_eval.h poly8_eval.h poly8_eval_low_mem.h
LIBS=-llapack
sp8cc.so: $(HEADER_FILES) cwrappers.cc
	$(CXX) -std=c++17 -fPIC -shared -o $@ cwrappers.cc $(LIBS)

clean:
	rm -f sp8cc.so