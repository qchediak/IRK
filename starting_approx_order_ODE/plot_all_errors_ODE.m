function result = plot_all_errors_ODE(stepsizes, errs_struct, params)
	% function result = plot_all_errors_ODE(stepsizes, errs_struct, params)
	%
	% Plots stepsizes vs. errors in log-log space and returns intercept
	% and slope of least-squares regression line.  Also prints
	% the least squares lines and the orders of the starting approximations.
	%
	% Inputs:
	%	stepsizes: Vector of stepsizes.
	%	errs_struct: Structured array of errors of starting
	%		approximations.
	%	params: Structured array with the following fields:
	%		problem_str: String describing the problem.
	%		method_str: String describing the method.
	%		starting_approx_str: Describes if the starting
	%			approximations are good or trivial.
	%		r: Ratio of stepsizes.
	%		s: Number of stages.
	%		iterative_scheme: Describes the iterative scheme
	%			used.
	%		error_plotting_Y: Error option for Y.  Either 'component' or
	%			'composite'.
	%
	% Output:
	%	result: Structured array of coefficients alpha, beta of the
	%		least-squares lines.

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	starting_approx_str = params.starting_approx_str;
	r = params.r;
	s = params.s;
	iterative_scheme = params.iterative_scheme;
	error_plotting_Y = params.error_plotting_Y;

	% unpack errs_struct
	Y0_errs = errs_struct.Y0_errs;

	% ----------------
	% Plot the errors
	% ----------------

	% For the pendulum problem, just write "pendulum", not "pendulum problem"
	if strcmpi(problem_str, 'pendulum')
		line2 =  sprintf('%s, %s, \n $s=%i$, $r=%.1f$', problem_str, method_str, s, r);
	else
		line2 =  sprintf('%s Problem, %s, \n $s=%i$, $r=%.1f$', problem_str, method_str, s, r);
	end

	fprintf('\n\n---------\nLeast Squares Lines\n---------\n')

	if strcmpi(error_plotting_Y, 'composite')

		line1_Y = '$\log_{10} \| Y_{n+1}^{(0)} - Y_{n+1} \|$ vs.~$\log_{10} h$';
		title_str_Y = sprintf('%s \n %s', line1_Y, line2);
		print_str_Y = 'Y0';

		plot_errs_params_Y.title_str = title_str_Y;
		plot_errs_params_Y.print_str = print_str_Y;

		[alpha_Y0, beta_Y0] = plot_initializer_errors(stepsizes, Y0_errs, plot_errs_params_Y);

	elseif strcmpi(error_plotting_Y, 'component') 

		for j=1:s
			line1_Y = sprintf('$\\log_{10} \\| Y_{n+1,%i}^{(0)} - Y_{n+1,%i} \\|$ vs.~$\\log_{10} h$', j, j);
			title_str_Y = sprintf('%s \n %s', line1_Y, line2);
			print_str_Y = sprintf('Y%i0',j);
			
			plot_errs_params_Y.title_str = title_str_Y;
			plot_errs_params_Y.print_str = print_str_Y;

			[alpha_Y0(j), beta_Y0(j)] = plot_initializer_errors(stepsizes, Y0_errs(:,j), plot_errs_params_Y);
		end

	end

	% make the result
	result.alpha_Y0 = alpha_Y0;
	result.beta_Y0 = beta_Y0;

	% -----------------------------------------------
	% Print the orders of the starting approximations
	% -----------------------------------------------

	fprintf('\n\n---------\nOrders of Starting Approximations\n---------\n')

	% Y Component
	if strcmpi(error_plotting_Y, 'composite')
		fprintf('\nOrder of starting approximation Y0: %0.1f', beta_Y0-1)
	elseif strcmpi(error_plotting_Y, 'component')
		for i=1:s
			fprintf('\nOrder of starting approximation Y%i0: %0.1f', i, beta_Y0(i)-1)
		end
	else
		error('Error (main): Invalid value of error_plotting_Y')
	end

	fprintf('\n\n')

end


