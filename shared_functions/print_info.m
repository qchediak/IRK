function print_info(method_str, problem_str, s, h, N, iter_struct, iterative_scheme_struct, number_nonlinear_systems, starting_approx_str, tol)
	% function print_info(method_str, problem_str, s, h, N, iter_struct, iterative_scheme_struct, number_nonlinear_systems, starting_approx_str, tol)
	%
	% Prints information on the method, problem, and iterations.
	%
	% Inputs:
	%	method_str: Describes the method.
	%	problem_str: Describes the problem.
	%	s: Number of stages.
	%	h: stepsize.
	%	N: Number of steps.
	%	iter_struct: Structured array of vector(s) of iterations.
	%	iterative_scheme_struct: Structured array of strings describing
	%		the iterative scheme(s) used.
	%	number_nonlinear_systems: The number of nonlinear systems being
	%		solved.
	%
	% Outputs:
	%	None.

	% unpack the structured arrays iter_struct and iterative_scheme_struct
	if number_nonlinear_systems==1
		iter_vec = iter_struct.iter_vec;
		iterative_scheme = iterative_scheme_struct.iterative_scheme;
	elseif number_nonlinear_systems==2
		iter1_vec = iter_struct.iter1_vec;
		iter2_vec = iter_struct.iter2_vec;
		iterative_scheme1 = iterative_scheme_struct.iterative_scheme1;
		iterative_scheme2 = iterative_scheme_struct.iterative_scheme2;
	end

	fprintf('\n--------------------------------------------')
	fprintf('\n%s, %s problem, s=%i, h=%01.e, N=%i', method_str, problem_str, s, h, N)
	fprintf('\n--------------------------------------------')

	if number_nonlinear_systems == 1

		iter_vec(1) = 0; % set to zero so the first step isn't counted in the sum

		iter_sum = sum(iter_vec);
		iter_ave = sum(iter_vec) / (N-1); % divide by N-1 since not counting the first step

		fprintf('\nStarting Approximations: %s', starting_approx_str)
		fprintf('\nTolerance: %0.1e', tol)
		fprintf('\nIterative Scheme: %s', iterative_scheme)
		fprintf('\nTotal Iterations (After First Step): %i', iter_sum);
		fprintf('\nAverage Iterations per step: %0.3f\n\n', iter_ave);

	elseif number_nonlinear_systems == 2

		iter1_vec(1) = 0; % set to zero so the first step isn't counted in the sum
		iter2_vec(1) = 0;

		iter1_sum = sum(iter1_vec);
		iter1_ave = sum(iter1_vec) / (N-1);

		iter2_sum = sum(iter2_vec);
		iter2_ave = sum(iter2_vec) / (N-1);

		iter12_sum = iter1_sum + iter2_sum;
		iter12_ave = iter12_sum / (2*(N-1));

		fprintf('\nStarting Approximations: %s', starting_approx_str)
		fprintf('\nTolerance: %0.1e\n', tol)

		fprintf('\nIterative Scheme (Nonlinear System 1): %s', iterative_scheme1)
		fprintf('\nTotal Iterations (After First Step) (Nonlinear System 1): %i', iter1_sum);
		fprintf('\nAverage Iterations per step (Nonlinear System 1): %0.3f', iter1_ave);

		fprintf('\n\nIterative Scheme (Nonlinear System 2): %s', iterative_scheme2)
		fprintf('\nTotal Iterations (After First Step) (Nonlinear System 2): %i', iter2_sum);
		fprintf('\nAverage Iterations per step (Nonlinear System 2): %0.3f\n', iter2_ave);

		fprintf('\nTotal Iterations (After First Step) (Nonlinear Systems 1 and 2): %i', iter12_sum);
		fprintf('\nAverage Iterations per step (Nonlinear Systems 1 and 2): %0.3f\n\n', iter12_ave);

	else

		error('Error (main): Invalid value of number_nonlinear_systems')

	end
end


