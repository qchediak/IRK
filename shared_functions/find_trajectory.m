function result = find_trajectory(iterative_scheme_struct, problem_str, method_str, starting_approx_str, h, s, tol, N)
	% function result = find_trajectory(iterative_scheme_struct, problem_str, method_str, starting_approx_str, h, s, tol, N)
	%
	% Solves an index 3 DAE and returns the values yns, zns, uns as a
	% structured array.
	%
	% Inputs:
	%	iterative_scheme_struct: Structured array of options for the
	%		iterative schemes used.  Has fields 'iterative_scheme',
	%		'iterative_scheme1', and 'iterative_scheme2'
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
	%	result: Structured array with fields 'yns', 'zns', 'uns', and
	%		'iterations_struct', the last of which gives information
	%		on the number of iterations used in solving the nonlinear
	%		systems.  Also includes the field 'E' for energy values if
	%		the problem is the pendulum.

	% ----------------------
	% Configure the method
	% ----------------------
	% Get values based on method_str and s
	method_result = get_values_method(method_str, s);

	% unpack method_result
	A_invertible = method_result.A_invertible;
	Ahat_invertible = method_result.Ahat_invertible;
	W = method_result.W;
	What = method_result.What;
	b = method_result.b;
	c = method_result.c;
	A = method_result.A;
	Ahat = method_result.Ahat;
	number_nonlinear_systems = method_result.number_nonlinear_systems;

	% unpack iterative_scheme_cell_array
	if number_nonlinear_systems==1
		iterative_scheme = iterative_scheme_struct.iterative_scheme;
	elseif number_nonlinear_systems==2
		iterative_scheme1 = iterative_scheme_struct.iterative_scheme1;
		iterative_scheme2 = iterative_scheme_struct.iterative_scheme2;
	else
		error('Error (find_trajectory): Invalid value of number_nonlinear_systems')
	end

	% ----------------------
	% Configure the problem
	% ----------------------

	problem_result = get_values_problem(problem_str);

	% unpack problem_result
	y0 = problem_result.y0;
	z0 = problem_result.z0;
	u0 = problem_result.u0;
	f = problem_result.f;
	k = problem_result.k;
	g = problem_result.g;
	dfdy = problem_result.dfdy;
	dfdz = problem_result.dfdz;
	dkdy = problem_result.dkdy;
	dkdz = problem_result.dkdz; 
	dkdu = problem_result.dkdu;
	dgdy = problem_result.dgdy;
	true_solution = problem_result.true_solution;

	y_true = true_solution.y_true;
	z_true = true_solution.z_true;
	u_true = true_solution.u_true;

	% Additional values for specific problems
	if strcmpi(problem_str, 'pendulum')
		g0 = problem_result.pendulum_struct.g0;
		l = problem_result.pendulum_struct.l;
		alpha0 = [0;g0];
	elseif strcmpi(problem_str, 'alternate pendulum')
		g0 = problem_result.pendulum_struct.g0;
		L = problem_result.pendulum_struct.L;
		m = problem_result.pendulum_struct.m;
		k0 = problem_result.pendulum_struct.k0;
		c0 = problem_result.pendulum_struct.c0;
		a0 = problem_result.pendulum_struct.a0;
		b0 = problem_result.pendulum_struct.b0;
		Minv = problem_result.pendulum_struct.Minv;
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		C = problem_result.C;
		dCdu = problem_result.dCdu;
		if strcmpi(problem_str, 'second nonlinear')
			M = problem_result.M;
		end
	end

	% Note that this assumes y0, z0, u0 are column vectors
	ny = size(y0,1); % dimensions of y
	nz = size(z0,1); % dimensions of z
	nu = size(u0,1); % dimensions of u

	% --------------------------
	% Iterative Scheme(s)
	% --------------------------

	% Check that the choice of iterative_scheme1 and iterative_scheme2 makes sense
	if number_nonlinear_systems==2
		if strcmpi(iterative_scheme1,'simplified newton 0') && ~strcmpi(method_str,'Lobatto IIIA-IIIB')
			error('Error (find_trajectory): The iterative scheme simplified newton 0 can only be used with the Lobatto IIIA-IIIB method')
		elseif (strcmpi(iterative_scheme1,'simplified newton 3') || strcmpi(iterative_scheme1, 'simplified newton index 3')) && (strcmpi(method_str,'Lobatto IIIA-IIIB') || strcmpi(method_str, 'projected Lobatto III'))
			error('Error (find_trajectory): The iterative schemes simplified newton 3 and simplified newton index 3 cannot be used with Lobatto IIIA-IIIB or projected Lobatto III.')
		elseif strcmpi(iterative_scheme2,'simplified newton 2') 
			error('Error (find_trajectory): The iterative scheme simplified newton 2 cannot be used with the second nonlinear system.')
		end
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
	SA_params.k = k;
	SA_params.g = g;
	SA_params.ny = ny;
	SA_params.nz = nz;
	SA_params.nu = nu;
	SA_params.method_str = method_str;
	SA_params.number_nonlinear_systems = number_nonlinear_systems;

	% Parameters to pass to the "solve" functions (functions which solve a
	% nonlinear system).
	% Note that params does NOT include h, N, y0, or z0.
	params.problem_str = problem_str;
	params.f = f;
	params.k = k;
	params.g = g;
	params.dfdy = dfdy;
	params.dfdz = dfdz;
	params.dkdy = dkdy;
	params.dkdz = dkdz;
	params.dkdu = dkdu;
	params.dgdy = dgdy;
	params.ny = ny;
	params.nz = nz;
	params.nu = nu;
	params.method_str = method_str;
	params.b = b;
	params.c = c;
	params.A = A;
	params.Ahat = Ahat;
	params.tol = tol;

	% Parameters to pass to the find_unp function.
	% Note that these do not change with the step number n.
	u_params.problem_str = problem_str;
	u_params.b = b;
	u_params.c = c;
	u_params.A_invertible = A_invertible;
	u_params.Ahat_invertible = Ahat_invertible;
	u_params.W = W;
	u_params.What = What;
	u_params.nu = nu;

	if strcmpi(problem_str, 'pendulum')
		u_params.g0 = g0;
		u_params.l = l;
	elseif strcmpi(problem_str, 'alternate pendulum')
		u_params.g0 = g0;
		u_params.L = L;
		u_params.m = m;
		u_params.k0 = k0;
		u_params.c0 = c0;
		u_params.a0 = a0;
		u_params.b0 = b0;
		u_params.Minv = Minv;
		u_params.dgdy = dgdy;
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		u_params.C = C;
		u_params.dCdu = dCdu;
		u_params.tol = tol;
	end

	% ---------
	% Main Part
	% ---------

	ynm = nan(ny, 1);
	znm = nan(nz, 1);
	unm = nan(nu, 1);

	Yn = nan(ny*s, 1);
	Zn = nan(nz*s, 1);
	Un = nan(nu*s, 1);

	yn = y0;
	zn = z0;
	un = u0;

	% matrices of yn and zn values
	yns = yn;
	zns = zn;
	uns = un;

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
			x0_struct = find_starting_approx(ynm, znm, yn, zn, un, Yn, Zn, Un, 'trivial', SA_params);
		elseif strcmpi(starting_approx_str, 'good')
			x0_struct = find_starting_approx(ynm, znm, yn, zn, un, Yn, Zn, Un, 'good', SA_params);
		else
			error('Error (find_trajectory): Invalid combination of n and starting_approx_str')
		end

		if number_nonlinear_systems == 1 % methods which have only one nonlinear system 

			% unpack x0_struct and make the starting approximation
			Ynp0 = x0_struct.Y0;
			Znp0 = x0_struct.Z0;
			Unp0 = x0_struct.U0;
			x0 = [Ynp0; Znp0; Unp0];

			[ynp znp Ynp Znp Unp iter] = solve_PRK_DAE_single_step(params, h, yn, zn, un, x0, iterative_scheme);

			iter_vec(n) = iter; % vector of iterations

		elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') 

			% unpack x0_struct and make the starting approximation
			Ynp0 = x0_struct.Y0;
			Znp0 = x0_struct.Z0;
			Unp0 = x0_struct.U0;
			znp0 = x0_struct.z0;
			Unps0 = x0_struct.Us0;

			x01 = [Ynp0; Znp0; Unp0];
			x02 = [znp0; Unps0];

			% First nonlinear system
			[Ybar Z Utilde iter1] = solve_first_system_lobattoIIIAIIIB(params, h, yn, zn, un, x01, iterative_scheme1);

			% Second nonlinear system
			[znp Us iter2] = solve_second_system_lobattoIIIAIIIB(params, h, yn, zn, Ybar, Z, Utilde, x02, iterative_scheme2);

			ynp = Ybar((s-2)*ny+1:(s-1)*ny); % ynp = Ys

			Ynp = [yn; Ybar];
			Znp = Z;
			Unp = [Utilde; Us];

			iter1_vec(n) = iter1;
			iter2_vec(n) = iter2;

		elseif strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')

			% unpack x0_struct and make the starting approximations
			Ynp0 = x0_struct.Y0;
			Znp0 = x0_struct.Z0;
			Unp0 = x0_struct.U0;
			znp0 = x0_struct.z0;
			Unpsp0 = x0_struct.Usp0;

			x01 = [Ynp0; Znp0; Unp0];
			x02 = [znp0; Unpsp0];

			% First nonlinear system.
			[ynp znp_not_used Ynp Znp Unp iter1] = solve_PRK_DAE_single_step(params, h, yn, zn, un, x01, iterative_scheme1);

			% Second nonlinear system
			[znp Unpsp iter2] = solve_second_system_projected(params, h, yn, zn, Ynp, Znp, Unp, x02, iterative_scheme2);

			iter1_vec(n) = iter1;
			iter2_vec(n) = iter2;

		elseif strcmpi(method_str, 'projected Lobatto III')

			% unpack x0_struct and make the starting approximations
			Yhat0 = x0_struct.Yhat0;
			Z0 = x0_struct.Z0;
			Uhat0 = x0_struct.Uhat0;
			znp0 = x0_struct.znp0;
			Us0 = x0_struct.Us0;

			x01 = [Yhat0; Z0; Uhat0];
			x02 = [znp0; Us0];

			% First nonlinear system.
			[Yhat Z Uhat iter1] = solve_first_system_lobattoIII(params, h, yn, zn, un, x01, iterative_scheme1);

			% Second nonlinear system
			[znp Us iter2] = solve_second_system_lobattoIII(params, h, yn, zn, Yhat, Z, Uhat, x02, iterative_scheme2);

			ynp = Yhat((s-1)*ny+1:s*ny);

			Ynp = [yn; Yhat(1:(s-1)*ny)]; % [Y1;...;Ys]
			Znp = Z;
			Unp = [Uhat(nu+1:s*nu); Us]; % [U1; ...; Us]

			iter1_vec(n) = iter1;
			iter2_vec(n) = iter2;

		else

			error('Error (find_trajectory): Invalid value of method_str')

		end

		% Choose a value for unp.  This is used in plotting and to set the
		% next (trivial) starting approximation(s).
		unp = find_unp(ynp, znp, Unp, un, u_params);

		% update the matrix of yns, zns
		yns = [yns ynp];
		zns = [zns znp];
		uns = [uns unp];

		if strcmpi(problem_str, 'pendulum')
			Enp = (1/2)*norm(znp)^2 + g0*(l+ynp(2)); % update energy
			E = [E Enp];
		end

		% update Yn, Zn, Un, ynm, znm, unm, yn, zn, un
		Yn = Ynp;
		Zn = Znp;
		Un = Unp; 

		ynm = yn;
		znm = zn;
		unm = un;

		yn = ynp;
		zn = znp;
		un = unp;

		% Debugging
		%if strcmpi(problem_str, 'second nonlinear')
		%	fprintf('\nn=%i, t=%0.2f', n, n*h)
		%	fprintf(' det(M) = %0.2f, cond(M)=%0.2f', det(M(ynp)), cond(M(ynp)))
		%end
	end	

	% -----------------------------------
	% Find the total number of iterations
	% -----------------------------------

	if number_nonlinear_systems==1
		iterations_struct.iter_vec = iter_vec;
	elseif number_nonlinear_systems==2
		iterations_struct.iter1_vec = iter1_vec;
		iterations_struct.iter2_vec = iter2_vec;
	else
		error('Error (find_trajectory): Invalid value of number_nonlinear_systems')
	end

	result.yns = yns;
	result.zns = zns;
	result.uns = uns;
	result.iterations_struct = iterations_struct;

	if strcmpi(problem_str, 'pendulum')
		result.E = E;
	end
end



