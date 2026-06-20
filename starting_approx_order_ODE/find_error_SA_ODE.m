function result = find_error_SA_ODE(h, params, option)
	% function result = find_error_SA_ODE(h, params, option)
	% Finds the error in the starting approximations 
	% when solving an ODE.
	%
	% Only two steps are taken: The first with stepsize h,
	% and the second with stepsize rh.
	% 
	% The "option" argument determines the way in which the internal stages
	% Yn are found:
	%
	% option='FF' (Forward, Forward): 
	%	Two forward steps are taken.  If there is a manifold on which the
	%	exact solution to the ODE lies, the "consistent" value will be
	%	ynm.
	%
	% option='BF' (Backward, Forward) (Recommended): 
	%	A backwards step with stepsize -h is taken
	% 	from yn from which we find Yn_star, 
	% 	from which Yn is found by the relationship 
	%	Y_{ni} = Y_{n,s+1-i}^*.  Then, a step of size rh is taken from
	%	yn.	If there is a manifold on which the exact solution to the
	%	ODE lies, the "consistent" value will be yn. This is the
	%	recommended setting.
	%
	% option='BFF' (Backward, Forward, Forward):
	% 	A backwards step with stepsize -h is taken from yn, from which
	% 	we find the update step ynm.  Then, two forward steps of
	% 	stepsizes h and rh are taken beginning at ynm.
	%
	% Options BF and BFF both theoretically should give the same results and
	% allow the middle step yn to stay the same across varying
	% stepsizes and also be "consistent" (i.e., lie on the correct
	% manifolds if they exist).  The 'BFF' option is conceptually a more
	% direct test of the theory, but it may introduce a small numerical
	% error in going from yn -> ynm -> yn', in which case the values yn
	% and yn' differ by a small amount.  For this reason, the 'BF'
	% option is recommended.  	
	%
	% Inputs:
	% 	h: stepsize
	%	params: struct format, listed below.
	%	option: see above.
	%
	% Outputs: 
	% 	result: A structure array with fields 'steps' and 'err_vecs',
	% 	where
	%
	%		result.steps is a structure array of the steps, where
	%
	%			result.steps.yns = [y_{n-1}, y_n, y_{n+1}]
	%
	%		result.err_vecs has the fields
	%
	%			Y0_err_vec

	% unpack params
	problem_str = params.problem_str;
	method_str = params.method_str;
	iterative_scheme = params.iterative_scheme;
	r = params.r;
	f = params.f;
	dfdy = params.dfdy;
	ny = params.ny;
	b = params.b;
	c = params.c;
	A = params.A;
	bstar = params.bstar;
	cstar = params.cstar;
	Astar = params.Astar;
	tol = params.tol;
	A_invertible = params.A_invertible;
	Astar_invertible = params.Astar_invertible;
	W = params.W;
	Wstar = params.Wstar;
	reverse_method_str = params.reverse_method_str;

	if strcmpi(problem_str, 'pendulum')
		g0 = params.pendulum_params.g0;
		l = params.pendulum_params.l;
		alpha0 = [0;g0];
	elseif strcmpi(problem_str, 'two body')
		a = params.two_body_params.a;
		alpha0 = params.two_body_params.alpha0;
		T = params.two_body_params.T;
		e = params.two_body_params.e;
	end

	% Interpret the "initial" values differently depending on the option
	if strcmpi(option, 'BF') || strcmpi(option, 'BFF')
		yn = params.y0;
	elseif strcmpi(option, 'FF')
		ynm = params.y0;
	else
		error('Error (find_error_SA_ODE): Invalid value of option.')
	end

	s = size(A,1);

	% Parameters for the "solve" functions with positive h.  
	% This is a subset of params. By design this does NOT contain yn.
	solve_params_fields = {'problem_str', 'f', 'dfdy', 'ny', 'method_str', 'b', 'c', 'A', 'tol'};
	for i=1:numel(solve_params_fields)
		solve_params.(solve_params_fields{i}) = params.(solve_params_fields{i});
	end

	% Parameters for find_starting_approx
	% A subset of params except for h 
	SA_params_fields = {'b', 'c', 'A', 'r', 'f', 'ny', 'method_str'};
	for i=1:numel(SA_params_fields)
		SA_params.(SA_params_fields{i}) = params.(SA_params_fields{i});
	end

	SA_params.h = h;

	% Parameters for the "solve" function with reverse method coefficients
	reverse_params_fields = {'problem_str', 'f', 'dfdy', 'ny', 'tol'};
	for i=1:numel(reverse_params_fields)
		reverse_params.(reverse_params_fields{i}) = params.(reverse_params_fields{i});
	end
	
	reverse_params.b = bstar;
	reverse_params.c = cstar;
	reverse_params.A = Astar;
	reverse_params.method_str = reverse_method_str;

	% Parameters for find_starting_approx with reverse method coefficients and -h
	% If using trivial starting approximations in step 0, this only requires
	% changing method_str to reverse_method_str.
	reverse_SA_params_fields = {'r', 'f', 'ny'};
	for i=1:numel(reverse_SA_params_fields)
		reverse_SA_params.(reverse_SA_params_fields{i}) = params.(reverse_SA_params_fields{i});
	end

	reverse_SA_params.h = -h;
	reverse_SA_params.b = bstar;
	reverse_SA_params.c = cstar;
	reverse_SA_params.A = Astar;
	reverse_SA_params.method_str = reverse_method_str;

	if strcmpi(option, 'BF') || strcmpi(option, 'BFF')
		% ------
		% Step 0
		% ------

		% Find the value ynm which comes before the value yn
		% Note this step uses -h.

		% Make trivial starting approximations for step 0.  Note this uses parameters of the reverse method.
		Yn0_star = find_starting_approx_ODE(nan(ny,1),yn,nan(ny*s,1),'trivial',reverse_SA_params);

		% Note this uses -h.
		[ynm Yn_star iter0] = solve_IRK_ODE_single_step(reverse_params, -h, yn, Yn0_star, iterative_scheme);
	end

	if strcmpi(option, 'BF')
		% Find the internal stages Yn, using Y_{ni} = Y_{n,s+1-i}^*
		E = flip(eye(s));
		Yn = kron(E,eye(ny))*Yn_star;

	elseif strcmpi(option, 'BFF') || strcmpi(option, 'FF')
		% ------
		% Step 1
		% ------

		% find the trivial starting approximation
		% Note that because of the input 'trivial', this only really uses ynm and s
		Y0_step1 = find_starting_approx_ODE(nan(ny,1), ynm, nan(ny*s,1), 'trivial', SA_params);

		[yn Yn iter] = solve_IRK_ODE_single_step(solve_params, h, ynm, Y0_step1, iterative_scheme);
	end

	% ------
	% Step 2
	% ------

	% find the good starting approximation
	Y0_step2 = find_starting_approx_ODE(ynm, yn, Yn, 'good', SA_params);

	[ynp Ynp iter] = solve_IRK_ODE_single_step(solve_params, r*h, yn, Y0_step2, iterative_scheme);

	% find the error
	Ynp0 = Y0_step2;
	Y0_err_vec = abs(Ynp0 - Ynp); 

	err_vecs.Y0_err_vec = Y0_err_vec;

	yns = [ynm yn ynp];
	steps.yns = yns;

	result.steps = steps;
	result.err_vecs = err_vecs;
end


