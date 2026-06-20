function result = find_error_SA(h, params, option)
	% function result = find_error_SA(h, params, option)
	% Finds the error in the starting approximations 
	% when solving an index 3 DAE.
	%
	% Only two steps are taken: The first with stepsize h,
	% and the second with stepsize rh.
	% 
	% The "option" argument determines the way in which the internal stages
	% Yn, Zn, Un are found:
	%
	% option='FF' (Forward, Forward): 
	%	Two forward steps are taken.  The consistent value will be
	%	(ynm,znm,unm).  This means that the step (yn,zn,un) will in
	%	general NOT be consistent.  The advantage to this option is it
	%	does not require implementing the reverse method.
	%
	% option='BF' (Backward, Forward) (Recommended): 
	%	A backwards step with stepsize -h is taken
	% 	from (yn,zn,un), from which we find Yn_star, Zn_star, Un_star,
	% 	from which Yn, Zn, Un are found by the relationship X_{ni} =
	% 	X_{n,s+1-i}^*.  Then, a step of size rh is taken from (yn,zn,un).
	%	The consistent value will be (yn,zn,un). This is the recommended setting.
	%
	% option='BFF' (Backward, Forward, Forward):
	% 	A backwards step with stepsize -h is taken from (yn,zn,un), from which
	% 	we find the update step (ynm,znm,unm).  Then, two forward steps of
	% 	stepsizes h and rh are taken beginning at (ynm,znm,unm).
	%
	% Options BF and BFF both theoretically should give the same results and
	% allow the middle step (yn,zn,un) to stay the same across varying
	% stepsizes and also be consistent (i.e., lie on the manifolds).  The
	% 'BFF' option is conceptually a more direct test of the theory, but it
	% may introduce a small numerical error in going from (yn,zn,un) ->
	% (ynm,znm,unm) -> (yn',zn',un'), in which case the values (yn,zn,un) and
	% (yn',zn',un') differ by a small amount.  For this reason, the 'BF'
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
	%			result.steps.zns = [z_{n-1}, z_n, z_{n+1}]
	%			result.steps.uns = [u_{n-1}, u_n, u_{n+1}].
	%
	%		result.err_vecs has three possibilities:
	%
	%		1. methods which solve only one nonlinear system:
	%			In this case result.err_vecs has the fields
	%			Y0_err_vec, Z0_err_vec, U0_err_vec
	%		2. Lobatto IIIA-IIIB: 
	%			In this case result.err_vecs has the fields
	%			Y0_err_vec, Z0_err_vec, U0_err_vec, znp0_err_vec,
	%			Unps0_err_vec
	%		3. Projected Lobatto IIIC: 
	%			In this case result.err_vecs has the fields
	%			Y0_err_vec, Z0_err_vec, U0_err_vec, znp0_err_vec,
	%			Unpsp0_err_vec
	%
	%		The possible components of err_vec are 
	%	
	%		Ynp0_err_vec: Vector of errors in the initializer Y_{n+1}^{(0)}
	%		Znp0_err_vec: Vector of errors in the initializer Z_{n+1}^{(0)}
	%		Unp0_err_vec: Vector of errors in the initializer U_{n+1}^{(0)}
	%		znp0_err_vec: Vector of errors in the initializer z_{n+1}^{(0)}
	%		Unps0_err_vec: Vector of errors in the initializer U_{n+1,s}^{(0)}
	%		Unpsp0_err_vec: Vector of errors in the initializer U_{n+1,s+1}^{(0)}

	% unpack params
	starting_approx_str = params.starting_approx_str;
	problem_str = params.problem_str;
	method_str = params.method_str;
	iterative_scheme_struct = params.iterative_scheme_struct;
	r = params.r;
	f = params.f;
	k = params.k;
	g = params.g;
	dfdy = params.dfdy;
	dfdz = params.dfdz;
	dkdy = params.dkdy;
	dkdz = params.dkdz;
	dkdu = params.dkdu;
	dgdy = params.dgdy;
	ny = params.ny;
	nz = params.nz;
	nu = params.nu;
	b = params.b;
	c = params.c;
	A = params.A;
	Ahat = params.Ahat;
	bstar = params.bstar;
	cstar = params.cstar;
	Astar = params.Astar;
	Ahatstar = params.Ahatstar;
	tol = params.tol;
	A_invertible = params.A_invertible;
	Ahat_invertible = params.Ahat_invertible;
	Astar_invertible = params.Astar_invertible;
	Ahatstar_invertible = params.Ahatstar_invertible;
	W = params.W;
	What = params.What;
	Wstar = params.Wstar;
	Whatstar = params.Whatstar;
	number_nonlinear_systems = params.number_nonlinear_systems;
	reverse_method_str = params.reverse_method_str;

	% Additional values for specific problems
	if strcmpi(problem_str, 'pendulum')
		g0 = params.pendulum_params.g0;
		l = params.pendulum_params.l;
		alpha0 = [0;g0];
	elseif strcmpi(problem_str, 'alternate pendulum')
		g0 = params.pendulum_params.g0;
		L = params.pendulum_params.L;
		m = params.pendulum_params.m;
		k0 = params.pendulum_params.k0;
		c0 = params.pendulum_params.c0;
		a0 = params.pendulum_params.a0;
		b0 = params.pendulum_params.b0;
		Minv = params.pendulum_params.Minv;
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		C = params.nonlinear_params.C;
		dCdu = params.nonlinear_params.dCdu;
		if strcmpi(problem_str, 'second nonlinear')
			M = params.nonlinear_params.M;
		end
	end

	% Interpret the "initial" values differently depending on the option
	if strcmpi(option, 'BF') || strcmpi(option, 'BFF')
		yn = params.y0;
		zn = params.z0;
		un = params.u0;
	elseif strcmpi(option, 'FF')
		ynm = params.y0;
		znm = params.z0;
		unm = params.u0;
	else
		error('Error (find_error_SA): Invalid value of option.')
	end

	% unpack iterative_scheme_struct
	if number_nonlinear_systems==1
		iterative_scheme = iterative_scheme_struct.iterative_scheme;
	elseif number_nonlinear_systems==2
		iterative_scheme1 = iterative_scheme_struct.iterative_scheme1;
		iterative_scheme2 = iterative_scheme_struct.iterative_scheme2;
	else
		error('Error (find_error_SA): Invalid value of number_nonlinear_systems')
	end

	s = size(A,1);

	% Parameters for the "solve" functions with positive h.  
	% This is a subset of params. By design this does NOT contain (yn,zn,un).
	solve_params_fields = {'problem_str', 'f', 'k', 'g', 'dfdy', 'dfdz', 'dkdy', 'dkdz', 'dkdu', 'dgdy', 'ny', 'nz', 'nu', 'method_str', 'b', 'c', 'A', 'Ahat', 'tol'};
	for i=1:numel(solve_params_fields)
		solve_params.(solve_params_fields{i}) = params.(solve_params_fields{i});
	end

	% Parameters for find_starting_approx
	% A subset of params except for h and number_nonlinear_systems.
	SA_params_fields = {'b', 'c', 'A', 'Ahat', 'r', 'f', 'k', 'ny', 'nz', 'nu', 'method_str'};
	for i=1:numel(SA_params_fields)
		SA_params.(SA_params_fields{i}) = params.(SA_params_fields{i});
	end

	SA_params.h = h;
	SA_params.number_nonlinear_systems = number_nonlinear_systems;

	% Parameters for the "solve" function with reverse method coefficients
	reverse_params_fields = {'problem_str', 'f', 'k', 'g', 'dfdy', 'dfdz', 'dkdy', 'dkdz', 'dkdu', 'dgdy', 'ny', 'nz', 'nu', 'tol'};
	for i=1:numel(reverse_params_fields)
		reverse_params.(reverse_params_fields{i}) = params.(reverse_params_fields{i});
	end
	
	reverse_params.b = bstar;
	reverse_params.c = cstar;
	reverse_params.A = Astar;
	reverse_params.Ahat = Ahatstar;
	reverse_params.method_str = reverse_method_str;

	% Parameters for find_starting_approx with reverse method coefficients and -h
	% If using trivial starting approximations in step 0, this only requires
	% changing method_str to reverse_method_str.
	reverse_SA_params_fields = {'r', 'f', 'k', 'ny', 'nz', 'nu'};
	for i=1:numel(reverse_SA_params_fields)
		reverse_SA_params.(reverse_SA_params_fields{i}) = params.(reverse_SA_params_fields{i});
	end

	reverse_SA_params.h = -h;
	reverse_SA_params.number_nonlinear_systems = number_nonlinear_systems;
	reverse_SA_params.b = bstar;
	reverse_SA_params.c = cstar;
	reverse_SA_params.A = Astar;
	reverse_SA_params.Ahat = Ahatstar;
	reverse_SA_params.method_str = reverse_method_str;

	% parameters for find_unp (forward direction)
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

	% parameters for find_unp (reverse direction)
	reverse_u_params.problem_str = problem_str;
	reverse_u_params.b = bstar;
	reverse_u_params.c = cstar;
	reverse_u_params.A_invertible = Astar_invertible;
	reverse_u_params.Ahat_invertible = Ahatstar_invertible;
	reverse_u_params.W = Wstar;
	reverse_u_params.What = Whatstar;
	reverse_u_params.nu = nu;

	if strcmpi(problem_str, 'pendulum')
		reverse_u_params.g0 = g0;
		reverse_u_params.l = l;
	elseif strcmpi(problem_str, 'alternate pendulum')
		reverse_u_params.g0 = g0;
		reverse_u_params.L = L;
		reverse_u_params.m = m;
		reverse_u_params.k0 = k0;
		reverse_u_params.c0 = c0;
		reverse_u_params.a0 = a0;
		reverse_u_params.b0 = b0;
		reverse_u_params.Minv = Minv;
		reverse_u_params.dgdy = dgdy;
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		u_params.C = C;
		u_params.dCdu = dCdu;
		u_params.tol = tol;
	end

	if strcmpi(option, 'BF') || strcmpi(option, 'BFF')
		% ------
		% Step 0
		% ------

		% Find the value (ynm,znm,unm) which comes before the consistent value (yn,zn,un)
		% Note this step uses -h.

		% Make trivial starting approximations for step 0.  Note this uses parameters of the reverse method.
		x0_struct_step0 = find_starting_approx(nan(ny,1),nan(nz,1),yn,zn,un,nan(ny*s,1),nan(nz*s,1),nan(nu*s,1),'trivial',reverse_SA_params);

		if number_nonlinear_systems==1

			% unpack x0_struct_step0 and make the starting approximation
			Yn0_star = x0_struct_step0.Y0;
			Zn0_star = x0_struct_step0.Z0;
			Un0_star = x0_struct_step0.U0;
			x0_step0 = [Yn0_star; Zn0_star; Un0_star];

			% Note this uses -h.
			[ynm znm Yn_star Zn_star Un_star iter0] = solve_PRK_DAE_single_step(reverse_params, -h, yn, zn, un, x0_step0, iterative_scheme);

			% set the value of Uns (may be used for un below)
			Uns_star = Un_star((s-1)*nu+1:s*nu);
		
		elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')

			% Unpack x0_struct_step0 and make the starting approximation.
			Yn0_star = x0_struct_step0.Y0;
			Zn0_star = x0_struct_step0.Z0;
			Un0_star = x0_struct_step0.U0;
			znm0 = x0_struct_step0.z0;
			Uns0_star = x0_struct_step0.Us0;

			x01_step0 = [Yn0_star; Zn0_star; Un0_star];
			x02_step0 = [znm0; Uns0_star];

			% First nonlinear system.  Uses -h
			[Ynbar_star Zn_star Untilde_star iter01] = solve_first_system_lobattoIIIAIIIB(reverse_params, -h, yn, zn, un, x01_step0, iterative_scheme1);

			% Second nonlinear system.  Uses -h
			[znm Uns_star iter02] = solve_second_system_lobattoIIIAIIIB(reverse_params, -h, yn, zn, Ynbar_star, Zn_star, Untilde_star, x02_step0, iterative_scheme2);

			ynm = Ynbar_star((s-2)*ny+1:(s-1)*ny); % ynm = Yns_star

			Yn_star = [yn; Ynbar_star];
			Un_star = [Untilde_star; Uns_star];

		elseif strcmpi(method_str, 'projected Lobatto IIIC')

			% unpack x0_struct_step0 and make the starting approximation
			Yhat0_star = x0_struct_step0.Yhat0;
			Zn0_star = x0_struct_step0.Z0;
			Uhat0_star = x0_struct_step0.Uhat0;
			znm0 = x0_struct_step0.znp0;
			Uns0_star = x0_struct_step0.Us0;

			x01_step0 = [Yhat0_star; Zn0_star; Uhat0_star];
			x02_step0 = [znm0; Uns0_star];

			% first_nonlinear_system of the reverse method, Lobatto III
			[Yhat_star Zn_star Uhat_star iter01] = solve_first_system_lobattoIII(reverse_params, -h, yn, zn, un, x01_step0, iterative_scheme1);

			% Second nonlinear system of the reverse method, Lobatto III
			[znm Uns_star iter02] = solve_second_system_lobattoIII(reverse_params, -h, yn, zn, Yhat_star, Zn_star, Uhat_star, x02_step0, iterative_scheme2);

			Yn_star = [yn; Yhat_star(1:(s-1)*ny)]; % [Y1_star; ...; Ys_star]
			Un_star = [Uhat_star(nu+1:s*nu); Uns_star]; % [U1_star; ...; Us_star]

			ynm = Yhat_star((s-1)*ny+1:s*ny);

		else

			error('Error (find_error_SA): Invalid value of method_str')

		end

		unm = find_unp(ynm,znm,Un_star,un,reverse_u_params);

	end

	if strcmpi(option, 'BF')
		% Find the internal stages Yn, Zn, Un using X_{ni} = X_{n,s+1-i}^*
		E = flip(eye(s));
		Yn = kron(E,eye(ny))*Yn_star;
		Zn = kron(E,eye(nz))*Zn_star;
		Un = kron(E,eye(nu))*Un_star;

	elseif strcmpi(option, 'BFF') || strcmpi(option, 'FF')
		% ------
		% Step 1
		% ------

		% find the trivial starting approximation
		% Note that because of the input 'trivial', this only really uses ynm, znm, unm, and s
		x0_struct_step1 = find_starting_approx(nan(ny,1), nan(nz,1), ynm, znm, unm, nan(ny*s,1), nan(nz*s,1), nan(nu*s,1), 'trivial', SA_params);

		if number_nonlinear_systems==1

			% unpack x0_struct_step1 and make the starting approximation
			Yn0 = x0_struct_step1.Y0;
			Zn0 = x0_struct_step1.Z0;
			Un0 = x0_struct_step1.U0;
			x0_step1 = [Yn0; Zn0; Un0];

			% first_nonlinear_system
			[yn zn Yn Zn Un iter] = solve_PRK_DAE_single_step(solve_params, h, ynm, znm, unm, x0_step1, iterative_scheme);

			% set the value of Uns (may be used for un below)
			Uns = Un((s-1)*nu+1:s*nu);
		
		elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')

			% Unpack x0_struct_step1 and make the starting approximation.
			Yn0 = x0_struct_step1.Y0;
			Zn0 = x0_struct_step1.Z0;
			Un0 = x0_struct_step1.U0;
			zn0 = x0_struct_step1.z0;
			Uns0 = x0_struct_step1.Us0;

			x01_step1 = [Yn0; Zn0; Un0];
			x02_step1 = [zn0; Uns0];

			% First nonlinear system
			[Ynbar Zn Untilde iter1] = solve_first_system_lobattoIIIAIIIB(solve_params, h, ynm, znm, unm, x01_step1, iterative_scheme1);

			% Second nonlinear system
			[zn Uns iter2] = solve_second_system_lobattoIIIAIIIB(solve_params, h, ynm, znm, Ynbar, Zn, Untilde, x02_step1, iterative_scheme2);

			yn = Ynbar((s-2)*ny+1:(s-1)*ny); % yn = Yns

			Yn = [ynm; Ynbar];
			Un = [Untilde; Uns];

		elseif strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')

			% unpack x0_struct_step1 and make the starting approximation
			Yn0 = x0_struct_step1.Y0;
			Zn0 = x0_struct_step1.Z0;
			Un0 = x0_struct_step1.U0;
			zn0 = x0_struct_step1.z0;
			Unsp0 = x0_struct_step1.Usp0;

			x01_step1 = [Yn0; Zn0; Un0];
			x02_step1 = [zn0; Unsp0];

			% first_nonlinear_system
			[yn zn_notused Yn Zn Un iter1] = solve_PRK_DAE_single_step(solve_params, h, ynm, znm, unm, x01_step1, iterative_scheme1);

			% Second nonlinear system
			[zn Unsp iter2] = solve_second_system_projected(solve_params, h, ynm, znm, Yn, Zn, Un, x02_step1, iterative_scheme2);

			% set the value of Uns (may be used for un below)
			Uns = Un((s-1)*nu+1:s*nu);

		else

			error('Error (find_error_SA): Invalid value of method_str')

		end

		un = find_unp(yn,zn,Un,unm,u_params);

	end

	% ------
	% Step 2
	% ------

	% Find the 'good' starting approximation or 'interpolation' starting approximation.
	% You could also use 'trivial' here, although there would be no point.
	x0_struct_step2 = find_starting_approx(ynm, znm, yn, zn, un, Yn, Zn, Un, starting_approx_str, SA_params);

	if number_nonlinear_systems==1

		% unpack x0_struct_step2 and make the starting approximation
		Ynp0 = x0_struct_step2.Y0;
		Znp0 = x0_struct_step2.Z0;
		Unp0 = x0_struct_step2.U0;
		x0_step2 = [Ynp0; Znp0; Unp0];

		% first_nonlinear_system
		[ynp znp Ynp Znp Unp iter] = solve_PRK_DAE_single_step(solve_params, r*h, yn, zn, un, x0_step2, iterative_scheme);

		% set the value of Unps (may be used for unp below)
		Unps = Unp(end-nu+1:end);

		% find the error
		Y0_err_vec = abs(Ynp0 - Ynp); 
		Z0_err_vec = abs(Znp0 - Znp);
		U0_err_vec = abs(Unp0 - Unp); 

		err_vecs.Y0_err_vec = Y0_err_vec;
		err_vecs.Z0_err_vec = Z0_err_vec;
		err_vecs.U0_err_vec = U0_err_vec;

	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')

		% unpack x0_struct_step2 and make the starting approximation
		Ynp0 = x0_struct_step2.Y0;
		Znp0 = x0_struct_step2.Z0;
		Unp0 = x0_struct_step2.U0;
		znp0 = x0_struct_step2.z0;
		Unps0 = x0_struct_step2.Us0;

		% Note that for Lobatto IIIA-IIIB Ynp0=(Y20,...,Ys0).
		x01_step2 = [Ynp0; Znp0; Unp0];
		x02_step2 = [znp0; Unps0];

		% First nonlinear system
		[Ynpbar Znp Unptilde iter3] = solve_first_system_lobattoIIIAIIIB(solve_params, r*h, yn, zn, un, x01_step2, iterative_scheme1);

		% Second nonlinear system
		[znp Unps iter4] = solve_second_system_lobattoIIIAIIIB(solve_params, r*h, yn, zn, Ynpbar, Znp, Unptilde, x02_step2, iterative_scheme2);

		ynp = Ynpbar((s-2)*ny+1:(s-1)*ny); % ynp = Ynps

		Ynp = [yn; Ynpbar];
		Unp = [Unptilde; Unps];

		% find the error
		Y0_err_vec = abs(Ynp0 - Ynpbar); % recall Ynp0 = [Y_2^(0); ...; Y_s^(0)]
		Y0_err_vec = [zeros(ny,1); Y0_err_vec]; % Set error of Y_1^0 equal to zero.  Keeps dimensions consistent.
		Z0_err_vec = abs(Znp0 - Znp);
		U0_err_vec = abs(Unp0 - Unptilde); % recall Unp0 = [U_1^(0); ...; U_{s-1}^(0)]

		znp0_err_vec = abs(znp0 - znp);
		Unps0_err_vec = abs(Unps0 - Unps);

		err_vecs.Y0_err_vec = Y0_err_vec;
		err_vecs.Z0_err_vec = Z0_err_vec;
		err_vecs.U0_err_vec = U0_err_vec;
		err_vecs.znp0_err_vec = znp0_err_vec;
		err_vecs.Unps0_err_vec = Unps0_err_vec;

	elseif strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')

		% unpack x0_struct_step2 and make the starting approximation
		Ynp0 = x0_struct_step2.Y0;
		Znp0 = x0_struct_step2.Z0;
		Unp0 = x0_struct_step2.U0;
		znp0 = x0_struct_step2.z0;
		Unpsp0 = x0_struct_step2.Usp0;

		x01_step2 = [Ynp0; Znp0; Unp0];
		x02_step2 = [znp0; Unpsp0];

		% First nonlinear system
		[ynp znp_not_used Ynp Znp Unp iter3] = solve_PRK_DAE_single_step(solve_params, r*h, yn, zn, un, x01_step2, iterative_scheme1);

		% Second nonlinear system
		[znp Unpsp iter4] = solve_second_system_projected(solve_params, r*h, yn, zn, Ynp, Znp, Unp, x02_step2, iterative_scheme2);

		% set the value of Unps (may be used for unp below)
		Unps = Unp((s-1)*nu+1:s*nu);

		% find the error
		Y0_err_vec = abs(Ynp0 - Ynp);
		Z0_err_vec = abs(Znp0 - Znp);
		U0_err_vec = abs(Unp0 - Unp);

		znp0_err_vec = abs(znp0 - znp);
		Unpsp0_err_vec = abs(Unpsp0 - Unpsp);

		err_vecs.Y0_err_vec = Y0_err_vec;
		err_vecs.Z0_err_vec = Z0_err_vec;
		err_vecs.U0_err_vec = U0_err_vec;
		err_vecs.znp0_err_vec = znp0_err_vec;
		err_vecs.Unpsp0_err_vec = Unpsp0_err_vec;

	else

		error('Error (find_error_SA): Invalid value of method_str')

	end

	unp = find_unp(ynp,znp,Unp,un,u_params);

	yns = [ynm yn ynp];
	zns = [znm zn znp];
	uns = [unm un unp];

	steps.yns = yns;
	steps.zns = zns;
	steps.uns = uns;

	result.steps = steps;
	result.err_vecs = err_vecs;
end

