function result = plot_all_errors(stepsizes, errs_struct, params)
	% function result = plot_all_errors(stepsizes, errs_struct, params)
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
	%		number_nonlinear_systems: The number of nonlinear systems to
	%			solve.
	%		r: Ratio of stepsizes.
	%		s: Number of stages.
	%		iterative_scheme_struct: Describes the iterative_schemes
	%			used.
	%		error_plotting_struct: Structured array describing the error
	%			plotting options.
	%
	% Output:
	%	result: Structured array of coefficients alpha, beta of the
	%		least-squares lines.

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	starting_approx_str = params.starting_approx_str;
	number_nonlinear_systems = params.number_nonlinear_systems;
	r = params.r;
	s = params.s;
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

	% unpack errs_struct
	Y0_errs = errs_struct.Y0_errs;
	Z0_errs = errs_struct.Z0_errs;
	U0_errs = errs_struct.U0_errs;
	if number_nonlinear_systems==2 && error_plotting_option==1
		znp0_errs = errs_struct.znp0_errs;
		if strcmpi(method_str,'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
			Unpsp0_errs = errs_struct.Unpsp0_errs;
		end
	end 

	% ----------------
	% Plot the errors
	% ----------------

	% For the pendulum problem, just write "pendulum", not "pendulum problem"
	% For the alternate pendulum problem, write "pendulum-spring" instead of "alternate pendulum problem"
	if strcmpi(problem_str, 'pendulum')
		line2 =  sprintf('%s, %s, \n $s=%i$, $r=%.1f$', problem_str, method_str, s, r);
	elseif strcmpi(problem_str, 'alternate pendulum')
		line2 =  sprintf('Pendulum-Spring, %s, \n $s=%i$, $r=%.1f$', method_str, s, r);
	else
		line2 =  sprintf('%s Problem, %s, \n $s=%i$, $r=%.1f$', problem_str, method_str, s, r);
	end
		

	fprintf('\n\n---------\nLeast Squares Lines\n---------\n')

	if error_plotting_option==1

		% Y component
		if strcmpi(error_plotting_Y, 'composite')

			line1_Y = '$\log_{10} \| Y_{n+1}^{(0)} - Y_{n+1} \|$ vs.~$\log_{10} h$';
			title_str_Y = sprintf('%s \n %s', line1_Y, line2);
			print_str_Y = 'Y0';

			plot_errs_params_Y.title_str = title_str_Y;
			plot_errs_params_Y.print_str = print_str_Y;

			[alpha_Y0, beta_Y0] = plot_initializer_errors(stepsizes, Y0_errs, plot_errs_params_Y);

		elseif strcmpi(error_plotting_Y, 'component') 

			if number_nonlinear_systems==1 || strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
				for j=1:s
					line1_Y = sprintf('$\\log_{10} \\| Y_{n+1,%i}^{(0)} - Y_{n+1,%i} \\|$ vs.~$\\log_{10} h$', j, j);
					title_str_Y = sprintf('%s \n %s', line1_Y, line2);
					print_str_Y = sprintf('Y%i0',j);
					
					plot_errs_params_Y.title_str = title_str_Y;
					plot_errs_params_Y.print_str = print_str_Y;

					[alpha_Y0(j), beta_Y0(j)] = plot_initializer_errors(stepsizes, Y0_errs(:,j), plot_errs_params_Y);
				end
			elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
				alpha_Y0(1) = NaN;
				beta_Y0(1) = NaN;

				for j=2:s
					line1_Y = sprintf('$\\log_{10} \\| Y_{n+1,%i}^{(0)} - Y_{n+1,%i} \\|$ vs.~$\\log_{10} h$', j, j);
					title_str_Y = sprintf('%s \n %s', line1_Y, line2);
					print_str_Y = sprintf('Y%i0',j);

					plot_errs_params_Y.title_str = title_str_Y;
					plot_errs_params_Y.print_str = print_str_Y;

					[alpha_Y0(j), beta_Y0(j)] = plot_initializer_errors(stepsizes, Y0_errs(:,j), plot_errs_params_Y);
				end
			else
				error('Error (plot_all_errors): Invalid value of number_nonlinear_systems or method_str')
			end

		else

			error('Error (plot_all_errors): Invalid value of error_plotting_Y')

		end

		% Z component
		if strcmpi(error_plotting_Z, 'composite')

			line1_Z = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$';
			title_str_Z = sprintf('%s \n %s', line1_Z, line2);
			print_str_Z = 'Z0';

			plot_errs_params_Z.title_str = title_str_Z;
			plot_errs_params_Z.print_str = print_str_Z;

			[alpha_Z0, beta_Z0] = plot_initializer_errors(stepsizes, Z0_errs, plot_errs_params_Z);

		elseif strcmpi(error_plotting_Z, 'component') 

			for j=1:s
				line1_Z = sprintf('$\\log_{10} \\| Z_{n+1,%i}^{(0)} - Z_{n+1,%i} \\|$ vs.~$\\log_{10} h$', j, j);
				title_str_Z = sprintf('%s \n %s', line1_Z, line2);
				print_str_Z = sprintf('Z%i0',j);

				plot_errs_params_Z.title_str = title_str_Z;
				plot_errs_params_Z.print_str = print_str_Z;

				[alpha_Z0(j), beta_Z0(j)] = plot_initializer_errors(stepsizes, Z0_errs(:,j), plot_errs_params_Z);
			end

		else

			error('Error (plot_all_errors): Invalid value of error_plotting_Z')

		end

		% U component
		if strcmpi(error_plotting_U, 'composite')

			line1_U = '$\log_{10} \| U_{n+1}^{(0)} - U_{n+1} \|$ vs.~$\log_{10} h$';
			title_str_U = sprintf('%s \n %s', line1_U, line2);
			print_str_U = 'U0';

			plot_errs_params_U.title_str = title_str_U;
			plot_errs_params_U.print_str = print_str_U;

			[alpha_U0, beta_U0] = plot_initializer_errors(stepsizes, U0_errs, plot_errs_params_U);

		elseif strcmpi(error_plotting_U, 'component') 

			for j=1:s
				line1_U = sprintf('$\\log_{10} \\| U_{n+1,%i}^{(0)} - U_{n+1,%i} \\|$ vs.~$\\log_{10} h$', j, j);
				title_str_U = sprintf('%s \n %s', line1_U, line2);
				print_str_U = sprintf('U%i0',j);

				plot_errs_params_U.title_str = title_str_U;
				plot_errs_params_U.print_str = print_str_U;

				[alpha_U0(j), beta_U0(j)] = plot_initializer_errors(stepsizes, U0_errs(:,j), plot_errs_params_U);
			end

		else

			error('Error (plot_all_errors): Invalid value of error_plotting_U')

		end

		if number_nonlinear_systems==2
			% znp
			line1_znp = '$\log_{10} \| z_{n+1}^{(0)} - z_{n+1} \|$ vs.~$\log_{10} h$';
			title_str_znp = sprintf('%s \n %s', line1_znp, line2);
			print_str_znp = 'znp0';

			plot_errs_params_znp.title_str = title_str_znp;
			plot_errs_params_znp.print_str = print_str_znp;

			[alpha_znp0, beta_znp0] = plot_initializer_errors(stepsizes, znp0_errs, plot_errs_params_znp);
		end

		if strcmpi(method_str,'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
			% Unpsp
			line1_Unpsp = '$\log_{10} \| U_{n+1,s+1}^{(0)} - U_{n+1,s+1} \|$ vs.~$\log_{10} h$';
			title_str_Unpsp = sprintf('%s \n %s', line1_Unpsp, line2);
			print_str_Unpsp = 'Unpsp0';

			plot_errs_params_Unpsp.title_str = title_str_Unpsp;
			plot_errs_params_Unpsp.print_str = print_str_Unpsp;

			[alpha_Unpsp0, beta_Unpsp0] = plot_initializer_errors(stepsizes, Unpsp0_errs, plot_errs_params_Unpsp);
		end

	elseif error_plotting_option==2 && s==3 && strcmpi(method_str,'projected Lobatto IIIC')

		% Y component
		line1_Y1 = '$\log_{10} \| Y_{n+1,1}^{(0)} - Y_{n+1,1} \|$ vs.~$\log_{10} h$';
		title_str_Y1 = sprintf('%s \n %s', line1_Y1, line2);
		print_str_Y1 = 'Y10';
		
		plot_errs_params_Y1.title_str = title_str_Y1;
		plot_errs_params_Y1.print_str = print_str_Y1;

		[alpha_Y0(1), beta_Y0(1)] = plot_initializer_errors(stepsizes, Y0_errs(:,1), plot_errs_params_Y1);

		line1_Y2 = '$\log_{10} \| Y_{n+1,2}^{(0)} - Y_{n+1,2} \|$ vs.~$\log_{10} h$';
		title_str_Y2 = sprintf('%s \n %s', line1_Y2, line2);
		print_str_Y2 = 'Y20';

		plot_errs_params_Y2.title_str = title_str_Y2;
		plot_errs_params_Y2.print_str = print_str_Y2;

		[alpha_Y0(2), beta_Y0(2)] = plot_initializer_errors(stepsizes, Y0_errs(:,2), plot_errs_params_Y2);

		line1_Y3 = '$\log_{10} \| Y_{n+1,3}^{(0)} - Y_{n+1,3} \|$ vs.~$\log_{10} h$';
		title_str_Y3 = sprintf('%s \n %s', line1_Y3, line2);
		print_str_Y3 = 'Y30';

		plot_errs_params_Y3.title_str = title_str_Y3;
		plot_errs_params_Y3.print_str = print_str_Y3;

		[alpha_Y0(3), beta_Y0(3)] = plot_initializer_errors(stepsizes, Y0_errs(:,3), plot_errs_params_Y3);

		% Z component
		line1_Z = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$';
		title_str_Z = sprintf('%s \n %s', line1_Z, line2);
		print_str_Z = 'Z0';

		plot_errs_params_Z.title_str = title_str_Z;
		plot_errs_params_Z.print_str = print_str_Z;

		[alpha_Z0, beta_Z0] = plot_initializer_errors(stepsizes, Z0_errs, plot_errs_params_Z);

		% U component
		line1_U = '$\log_{10} \| U_{n+1}^{(0)} - U_{n+1} \|$ vs.~$\log_{10} h$';
		title_str_U = sprintf('%s \n %s', line1_U, line2);
		print_str_U = 'U0';

		plot_errs_params_U.title_str = title_str_U;
		plot_errs_params_U.print_str = print_str_U;

		[alpha_U0, beta_U0] = plot_initializer_errors(stepsizes, U0_errs, plot_errs_params_U);

	else

		error('Error (plot_all_errors): Invalid combination of error_plotting_option, s, and method_str.')

	end

	% make the result
	result.alpha_Y0 = alpha_Y0;
	result.alpha_Z0 = alpha_Z0;
	result.alpha_U0 = alpha_U0;

	result.beta_Y0 = beta_Y0;
	result.beta_Z0 = beta_Z0;
	result.beta_U0 = beta_U0;

	if number_nonlinear_systems==2 && error_plotting_option==1
		result.alpha_znp0 = alpha_znp0;
		result.beta_znp0 = beta_znp0;
		if strcmpi(method_str,'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
			result.alpha_Unpsp0 = alpha_Unpsp0;
			result.beta_Unpsp0 = beta_Unpsp0;
		end
	end

	% -----------------------------------------------
	% Print the orders of the starting approximations
	% -----------------------------------------------

	fprintf('\n\n---------\nOrders of Starting Approximations\n---------\n')

	if error_plotting_option==1
		% Y Component
		if strcmpi(error_plotting_Y, 'composite')
			fprintf('\nOrder of starting approximation Y0: %0.1f', beta_Y0-1)
		elseif strcmpi(error_plotting_Y, 'component')
			if number_nonlinear_systems==1 || strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
				for i=1:s
					fprintf('\nOrder of starting approximation Y%i0: %0.1f', i, beta_Y0(i)-1)
				end
			elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
				for i=2:s
					fprintf('\nOrder of starting approximation Y%i0: %0.1f', i, beta_Y0(i)-1)
				end
			else
				error('Error (main): Invalid value of method_str')
			end
		else
			error('Error (main): Invalid value of error_plotting_Y')
		end

		% Z Component
		if strcmpi(error_plotting_Z, 'composite')
			fprintf('\nOrder of starting approximation Z0: %0.1f', beta_Z0-1)
		elseif strcmpi(error_plotting_Z, 'component')
			for i=1:s
				fprintf('\nOrder of starting approximation Z%i0: %0.1f', i, beta_Z0(i)-1)
			end
		else
			error('Error (main): Invalid value of error_plotting_Z')
		end

		% U Component
		if strcmpi(error_plotting_U, 'composite')
			fprintf('\nOrder of starting approximation U0: %0.1f', beta_U0-1)
		elseif strcmpi(error_plotting_U, 'component')
			for i=1:s
				fprintf('\nOrder of starting approximation U%i0: %0.1f', i, beta_U0(i)-1)
			end
		else
			error('Error (main): Invalid value of error_plotting_U')
		end

		% znp
		if number_nonlinear_systems==2
			fprintf('\nOrder of starting approximation znp0: %0.1f', beta_znp0-1)
		end

		% Unpsp
		if strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
			fprintf('\nOrder of starting approximation Unpsp0: %0.1f', beta_Unpsp0-1)
		end
	elseif error_plotting_option==2
		% Y Component

		for i=1:s
			fprintf('\nOrder of starting approximation Y%i0: %0.1f', i, beta_Y0(i)-1)
		end

		% Z Component
		fprintf('\nOrder of starting approiximation Z0=(Z10,Z20,Z30,znp0): %0.1f', beta_Z0-1)

		% U Component
		fprintf('\nOrder of starting approiximation U0=(U10,U20,U30,U400): %0.1f', beta_U0-1)
	else
		error('Error (main): Invalid value of error_plotting_option')
	end

	fprintf('\n\n')

end

