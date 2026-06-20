function [znp Unps iter] = solve_second_system_lobattoIII(params, h, yn, zn, Yhat, Z, Uhat, x0, iterative_scheme)
	% function [znp Unps iter] = solve_second_system_lobattoIII(params, h, yn, zn, Yhat, Z, Uhat, x0, iterative_scheme)
	%
	% Solves the second nonlinear system in the projected Lobatto III
	% method.
	%
	% Inputs:
	%	params: Structured array consisting of parameters 
	%	h: step size.
	%	yn, zn: y- and z-components of the update step from the most
	%		recent step.
	%	Yhat, Z, Uhat: Vectors of internal stages already found in the
	%		first nonlinear system.  Here
	%		Yhat=(Y2,...,Ys,ynp), Z=(Z1,...,Zs), and 
	%		Uhat=(U0, U1,...,U_{s-1}).
	%	x0 = (znp0, Unps0): The starting approximation for the iterative
	%		solver for the second system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton': Uses the matrix DF(x0)
	%
	% Outputs:
	%	znp: The z-component of the update step
	% 	Unps: The stage U_{n+1,s}
	% 	iter: Number of iterations used by the iterative scheme.

	% unpack params
	problem_str = params.problem_str;
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
	method_str = params.method_str;
	b = params.b;
	c = params.c;
	A = params.A;
	Ahat = params.Ahat;
	tol = params.tol;

	s = size(A,1);

	% ---------------------------
	% Solve the nonlinear system
	% ---------------------------

	% The function T needs to be redefined at every step using Tparams.
	% The argument of T is x = [znp; Unpsp].

	% subset of params except for h, yn, zn, Yhat, Z, Uhat
	Tparams_fields = {'ny', 'nz', 'nu', 'b', 'f', 'k', 'dfdz', 'dgdy', 'dkdu'};
	for i=1:numel(Tparams_fields)
		Tparams.(Tparams_fields{i}) = params.(Tparams_fields{i});
	end

	Tparams.h = h;
	Tparams.yn = yn;
	Tparams.zn = zn;
	Tparams.Yhat = Yhat;
	Tparams.Z = Z;
	Tparams.Uhat = Uhat;

	T = @(x) T1_lobattoIII(x,Tparams);
	DT = @(x) DT1_lobattoIII(x,Tparams);

	if strcmpi(iterative_scheme, 'newton')
		[root iter] = newton(T,DT,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton') || strcmpi(iterative_scheme, 'simplified newton 1')
		[root iter] = simplified_newton(T,DT(x0),x0,tol);
	elseif strcmpi(iterative_scheme, 'kelley')

		combined = @(x) combined1(x,T,DT);
		tol_nsold = [1.d-10,1.d-10];
		params_nsold = [40, -1, 0.5, 0];
		%params_nsold = [40, -1, 0.5, 1];
		[root, it_hist, ierr, x_hist] = nsold(x0,combined,tol_nsold,params_nsold);
		iter = size(it_hist,1);
		
	else
		error('Error (solve_second_system_lobattoIII): Unknown iterative scheme');
	end

	% unpack root
	znp = root(1:nz);
	Unps = root(nz+1:nz+nu);
end


