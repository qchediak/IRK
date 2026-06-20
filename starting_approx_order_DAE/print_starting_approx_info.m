function print_starting_approx_info(stepsizes, params)
	% function print_starting_approx_info(stepsizes, params)
	%
	% Prints information relevant to finding the order of the starting
	% approximations for index 3 DAEs.

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	starting_approx_str = params.starting_approx_str;
	number_nonlinear_systems = params.number_nonlinear_systems;
	r = params.r;
	s = params.s;
	tol = params.tol;
	iterative_scheme_struct = params.iterative_scheme_struct;
	error_plotting_struct = params.error_plotting_struct;

	% unpack iterative_scheme_struct
	if number_nonlinear_systems==1
		iterative_scheme = iterative_scheme_struct.iterative_scheme;
	elseif number_nonlinear_systems==2
		iterative_scheme1 = iterative_scheme_struct.iterative_scheme1;
		iterative_scheme2 = iterative_scheme_struct.iterative_scheme2;
	else
		error('Error (print_starting_approx_info): Invalid value of number_nonlinear_systems.')
	end

	% unpack error_plotting_struct
	error_plotting_option = error_plotting_struct.error_plotting_option;
	error_plotting_Y = error_plotting_struct.error_plotting_Y;
	error_plotting_Z = error_plotting_struct.error_plotting_Z;
	error_plotting_U = error_plotting_struct.error_plotting_U;

	% Print some info
	fprintf('\n-------------')
	fprintf('\n%s Problem, %s, s=%i, r=%0.2f', problem_str, method_str, s, r)
	fprintf('\n-------------')
	fprintf('\nStarting Approximations: %s', starting_approx_str)
	fprintf('\nTolerance: %0.1e', tol)
	if number_nonlinear_systems==1
		fprintf('\nIterative Scheme: %s', iterative_scheme)
	elseif number_nonlinear_systems==2
		fprintf('\nIterative Scheme (Nonlinear System 1): %s', iterative_scheme1)
		fprintf('\nIterative Scheme (Nonlinear System 2): %s', iterative_scheme2)
	end
	fprintf('\nError Plotting Option: %i', error_plotting_option)
	if error_plotting_option==1
		fprintf('\nY, Z, U Plotting Options: %s, %s, %s', error_plotting_Y, error_plotting_Z, error_plotting_U)
	end
	fprintf('\nlog10(min(h)), log10(max(h)): %0.1f, %0.1f', log10(min(stepsizes)), log10(max(stepsizes)))

end

