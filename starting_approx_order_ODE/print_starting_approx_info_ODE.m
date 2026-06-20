function print_starting_approx_info_ODE(stepsizes, params)
	% function print_starting_approx_info_ODE(stepsizes, params)
	%
	% Prints information relevant to finding the order of the starting
	% approximations for ODEs.

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	starting_approx_str = params.starting_approx_str;
	r = params.r;
	s = params.s;
	tol = params.tol;
	iterative_scheme = params.iterative_scheme;
	error_plotting_Y = params.error_plotting_Y;

	% Print some info
	fprintf('\n-------------')
	fprintf('\n%s Problem, %s, s=%i, r=%0.2f', problem_str, method_str, s, r)
	fprintf('\n-------------')
	fprintf('\nStarting Approximations: %s', starting_approx_str)
	fprintf('\nTolerance: %0.1e', tol)
	fprintf('\nIterative Scheme: %s', iterative_scheme)
	fprintf('\nError Plotting Option: %s', error_plotting_Y)
	fprintf('\nlog10(min(h)), log10(max(h)): %0.1f, %0.1f', log10(min(stepsizes)), log10(max(stepsizes)))
end


