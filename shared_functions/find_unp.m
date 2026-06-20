function unp = find_unp(ynp, znp, Unp, un, params)
	% function unp = find_unp(ynp, znp, Unp, un, params)
	%
	% Finds the value update step for the u-component when solving an
	% index 3 DAE.
	%
	% Inputs:
	%	ynp, znp: Update steps for y- and z-components.
	%	Unp: Vector of internal stages for the u-component.
	% 	un: Previous step of the u-component.
	%	params: Structured array with the fields listed below.
	%		Note that nothing in params changes with the step number n.
	%		Fields of params are:
	%			problem_str: String describing the problem
	%			b: Column vector of weights.
	%			c: Column vector of nodes.
	%			A_invertible, Ahat_invertible: Boolean values describing
	%				invertibility of A and Ahat.
	%			W, What: Inverses of A and Ahat if they are nonsingular.
	%				If singular, then this is set to NaN.
	%			nu: Dimension of the u-component.
	%		For the 'pendulum'problem:
	%			g0: Gravitational constant 
	%			l: length of rod
	%		For the 'alternate pendulum' problem:
	%			g0: Gravitational constant 
	%			L: length of rod 
	%			m: mass
	%			k0: spring stiffness
	%			c0: damping coefficient
	%			a0, b0: Parts of the function k.
	%				k(y,z,u) = Minv*(a0(y,z) + b0(y)*u)
	%			Minv: Inverse of the mass matrix.
	%
	% Outputs:
	%	unp: Update step for the u-component.

	% unpack params
	problem_str = params.problem_str;
	b = params.b;
	c = params.c;
	A_invertible = params.A_invertible;
	Ahat_invertible = params.Ahat_invertible;
	W = params.W;
	What = params.What;
	nu = params.nu;

	if strcmpi(problem_str, 'pendulum')
		g0 = params.g0;
		l = params.l;
		alpha0 = [0;g0];
	elseif strcmpi(problem_str, 'alternate pendulum')
		g0 = params.g0;
		L = params.L;
		m = params.m;
		k0 = params.k0;
		c0 = params.c0;
		a0 = params.a0;
		b0 = params.b0;
		Minv = params.Minv;
		dgdy = params.dgdy;
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		C1 = params.C;
		dC1du = params.dCdu;
		tol = params.tol;
	end

	s = length(c); % number of stages

	% for initial values, Unp is sometimes passed as NaN
	if ~isnan(Unp)
		Us = Unp((s-1)*nu+1:s*nu); % U_{n+1,s}
	end

	% Choose a value for unp.  
	if strcmpi(problem_str, 'pendulum')
		% This is the value of unp which satisfies the hidden constraint.
		unp = (1/l^2) * (dot(znp,znp) - dot(ynp,alpha0)); 
	elseif strcmpi(problem_str, 'alternate pendulum')
		% This is the value of unp which satisfies the hidden constraint.
		LHS = @(y,z) dgdy(y)*Minv*b0(y);
		RHS = @(y,z) -L*z(3)^2*[cos(y(3)); sin(y(3))] - dgdy(y)*Minv*a0(y,z);
		unp = LHS(ynp,znp) \ RHS(ynp,znp);
	elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		% This is the value of unp which satisfies the hidden constraint.
		C = @(u) C1(ynp,znp,u);
		dC = @(u) dC1du(ynp,znp,u);
		[unp iter] = newton(C,dC,un,tol);
	elseif c(s)==1
		% Set unp = Us.  
		unp = Us;
	elseif A_invertible
		% find the value R(\infty)
		R = 1-b'*W*ones(s,1);
		unp = R*un + kron(b'*W,eye(nu))*Unp;
	elseif Ahat_invertible
		% find the value \widehat{R}(\infty)
		Rhat = 1-b'*What*ones(s,1);
		unp = Rhat*un + kron(b'*What,eye(nu))*Unp;
	else
		% At this point the best thing would be to solve the hidden
		% constraint analytically for unp.  The next best thing is
		% to set unp = Us.  No preset methods in this code should
		% enter this branch.
		unp = Us;
	end

end

