function result = find_trajectory_pODE(iterative_scheme, problem_str, method_str, starting_approx_str, h, s, tol, N)
	%function result = find_trajectory_pODE(iterative_scheme, problem_str, method_str, starting_approx_str, h, s, tol, N)
	%
	% Solves a partitioned ODE by a PRK method (b,c,A,Ahat) and returns
	% the values yns, zns, iterations, and possibly E (energy) as a
	% structured array.
	%
	% Inputs:
	%	iterative_scheme: String describing the iterative scheme.
	%	problem_str: String describing the problem.
	%	method_str: String describing the numerical method.
	%	starting_approx_str: Indicates whether 'good' or 'trivial'
	%		starting approximations should be used.
	%	h: Stepsize.  Constant steps are used.
	%	s: number of stages in the method
	%	tol: tolerance for root finder
	%	N: Number of steps to take.
	%
	% Outputs:
	%	result: Structured array with fields 'yns', 'zns', 'iter_vec',
	%	and possibly 'E'.

	% ----------------------
	% Configure the method
	% ----------------------
	% Get values based on method_str and s
	method_result = get_values_method(method_str, s);

	% unpack method_result
	b = method_result.b;
	c = method_result.c;
	A = method_result.A;
	Ahat = method_result.Ahat;

	% ----------------------
	% Configure the problem
	% ----------------------

	problem_result = get_values_problem_pODE(problem_str);

	% unpack problem_result
	y0 = problem_result.y0;
	z0 = problem_result.z0;
	f = problem_result.f;
	g = problem_result.g;
	dfdy = problem_result.dfdy;
	dfdz = problem_result.dfdz;
	dgdy = problem_result.dgdy;
	dgdz = problem_result.dgdz;

	% Additional values for specific problems
	if strcmpi(problem_str, 'pendulum')
		g0 = problem_result.pendulum_struct.g0;
		l = problem_result.pendulum_struct.l;
	elseif strcmpi(problem_str, 'two body')
		a = problem_result.two_body_struct.a;
		alpha0 = problem_result.two_body_struct.alpha0;
		T = problem_result.two_body_struct.T;
	end

	% Note that this assumes y0 and z0 are column vectors
	ny = size(y0,1); % dimensions of y
	nz = size(z0,1); % dimensions of z

	% --------------------------
	% Iterative Scheme(s)
	% --------------------------

	% Check that the choice of iterative_scheme1 and iterative_scheme2 makes sense
	if ~strcmpi(iterative_scheme, 'newton') && ~strcmpi(iterative_scheme, 'simplified newton') && ~strcmpi(iterative_scheme, 'simplified newton 1') && ~strcmpi(iterative_scheme, 'kelley')
		error('Error (find_trajectory_pODE): Invalid value of iterative_scheme.')
	end

	% -------------
	% Miscellaneous
	% -------------

	% Parameters needed to pass to find_starting_approx.
	% Note that this sets r=1 (constant steps)
	SA_params.h = h;
	SA_params.b = b;
	SA_params.c = c;
	SA_params.A = A;
	SA_params.Ahat = Ahat;
	SA_params.r = 1;
	SA_params.f = f;
	SA_params.g = g;
	SA_params.ny = ny;
	SA_params.nz = nz;
	SA_params.method_str = method_str;

	% Parameters to pass to the "solve" functions 
	% (functions which solve a nonlinear system).
	% Note that params does NOT include h, N, y0, or z0.
	params.problem_str = problem_str;
	params.f = f;
	params.g = g;
	params.dfdy = dfdy;
	params.dfdz = dfdz;
	params.dgdy = dgdy;
	params.dgdz = dgdz;
	params.ny = ny;
	params.nz = nz;
	params.method_str = method_str;
	params.b = b;
	params.c = c;
	params.A = A;
	params.Ahat = Ahat;
	params.tol = tol;

	% ---------
	% Main Part
	% ---------

	ynm = nan(ny, 1);
	znm = nan(nz, 1);
	Yn = nan(ny*s, 1);
	Zn = nan(nz*s, 1);
	yn = y0;
	zn = z0;

	% matrices of yn values
	yns = yn;
	zns = zn;

	if strcmpi(problem_str, 'pendulum')
		% initial energy, where potential energy is zero when the pendulum
		% is at its lowest point.
		E = (1/2)*norm(z0)^2 + g0*(l+y0(2)); 
	end

	% ----------
	% Main loop
	% ----------

	for n=1:N
		% Find the starting approximations.  Must be trivial in the first step.
		if strcmpi(starting_approx_str, 'trivial') || n==1
			[Y0 Z0] = find_starting_approx_pODE(ynm, znm, yn, zn, Yn, Zn, 'trivial', SA_params);
			x0 = [Y0; Z0];
		elseif strcmpi(starting_approx_str, 'good')
			[Y0 Z0] = find_starting_approx_pODE(ynm, znm, yn, zn, Yn, Zn, 'good', SA_params);
			x0 = [Y0; Z0];
		else
			error('Error (find_trajectory_pODE): Invalid combination of n and starting_approx_str')
		end

		[ynp znp Ynp Znp iter] = solve_PRK_pODE_single_step(params, h, yn, zn, x0, iterative_scheme);

		iter_vec(n) = iter; % vector of iterations

		% update the matrix of yns, zns
		yns = [yns ynp];
		zns = [zns znp];

		if strcmpi(problem_str, 'pendulum')
			Enp = (1/2)*norm(znp)^2 + g0*(l+ynp(2)); % update energy
			E = [E Enp];
		end

		% update Yn, Zn, ynm, znm, yn, and zn
		Yn = Ynp;
		Zn = Znp;
		ynm = yn;
		znm = zn;
		yn = ynp;
		zn = znp;
	end	

	result.yns = yns;
	result.zns = zns;
	result.iter_vec = iter_vec;

	if strcmpi(problem_str, 'pendulum')
		result.E = E;
	end
end



