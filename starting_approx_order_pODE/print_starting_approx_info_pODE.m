function print_starting_approx_info_pODE(stepsizes, params)
	% function print_starting_approx_info_pODE(stepsizes, params)
	%
	% Prints information relevant to finding the order of the starting
	% approximations for partitioned ODEs.

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	starting_approx_str = params.starting_approx_str;
	r = params.r;
	s = params.s;
	tol = params.tol;
	iterative_scheme = params.iterative_scheme;
	error_plotting_struct = params.error_plotting_struct;

	% unpack error_plotting_struct
	error_plotting_Y = error_plotting_struct.error_plotting_Y;
	error_plotting_Z = error_plotting_struct.error_plotting_Z;

	% Print some info
	fprintf('\n-------------')
	fprintf('\n%s Problem, %s, s=%i, r=%0.2f', problem_str, method_str, s, r)
	fprintf('\n-------------')
	fprintf('\nStarting Approximations: %s', starting_approx_str)
	fprintf('\nTolerance: %0.1e', tol)
	fprintf('\nIterative Scheme: %s', iterative_scheme)
	fprintf('\nY and Z Plotting Options: %s, %s', error_plotting_Y, error_plotting_Z)
	fprintf('\nlog10(min(h)), log10(max(h)): %0.1f, %0.1f', log10(min(stepsizes)), log10(max(stepsizes)))
end

