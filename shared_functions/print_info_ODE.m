function print_info_ODE(method_str, problem_str, s, h, N, iter_vec, iterative_scheme, starting_approx_str, tol)
	% function print_info_ODE(method_str, problem_str, s, h, N, iter_vec, iterative_scheme, starting_approx_str, tol)
	%
	% Prints information on the method, problem, and iterations.
	%
	% Inputs:
	%	method_str: Describes the method.
	%	problem_str: Describes the problem.
	%	s: Number of stages.
	%	h: stepsize.
	%	N: Number of steps.
	%	iter_vec: vector of iterations.
	%	iterative_scheme: String describing the iterative scheme used.
	%
	% Outputs:
	%	None.

	fprintf('\n--------------------------------------------')
	fprintf('\n%s, %s problem, s=%i, h=%01.e, N=%i', method_str, problem_str, s, h, N)
	fprintf('\n--------------------------------------------')

	iter_vec(1) = 0; % set to zero so the first step isn't counted in the sum

	iter_sum = sum(iter_vec);
	iter_ave = sum(iter_vec) / (N-1); % divide by N-1 since not counting the first step

	fprintf('\nStarting Approximations: %s', starting_approx_str)
	fprintf('\nTolerance: %0.1e', tol)
	fprintf('\nIterative Scheme: %s', iterative_scheme)
	fprintf('\nTotal Iterations (After First Step): %i', iter_sum);
	fprintf('\nAverage Iterations per step: %0.3f\n\n', iter_ave);
end



