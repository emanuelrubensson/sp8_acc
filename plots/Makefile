TOPDIR=..
plots:
	export PYTHONPATH="$(TOPDIR):$(PYTHONPATH)" && ./plot_sp8_acc.py
	export PYTHONPATH="$(TOPDIR):$(PYTHONPATH)" && ./plot_slope_at_mu.py
	export PYTHONPATH="$(TOPDIR):$(PYTHONPATH)" && ./plot_max_slope_at_mu.py
	export PYTHONPATH="$(TOPDIR):$(PYTHONPATH)" && ./plot_convergence.py

plot_convergence.pdf: plot_convergence.py
	export PYTHONPATH="$(TOPDIR):$(PYTHONPATH)" && ./plot_convergence.py
