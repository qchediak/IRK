function [znp Us iter] = solve_second_system_lobattoIIIAIIIB(params, h, yn, zn, Ybar, Z, Utilde, x0, iterative_scheme)
	% function [znp Us iter2] = solve_second_system_lobattoIIIAIIIB(params, h, yn, zn, Ybar, Z, Utilde, x0, iterative_scheme)
	%
	% Solve the second nonlinear system of the lobatto IIIA-IIIB method.
	%
	% Inputs:
	%	params: Structured array of parameters	.
	%	h: step size.
	%	(yn, zn): The most recent step. 	
	%	Ybar = (Y2,...Ys): Internal stages of the y-component
	%	Z = (Z1,...,Zs): Internal stages of the z-component
	%	Utilde = (U1,...,U_{s-1}): Internal stages of the u-component
	%	x0: initial guess.
	%	iterative_scheme: The iterative scheme used to solve the
	%		nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton': Uses the matrix DF(x0)
	%
	% Outputs:
	%	znp: The update step z_{n+1}
	%	Us: The final internal stage of the u-component
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

	% subset of params except for h, yn, zn, Ybar, Z, Utilde
	Tparams_fields = {'ny', 'nz', 'nu', 'f', 'k', 'g', 'dfdy', 'dfdz', 'dkdy', 'dkdz', 'dkdu', 'dgdy', 'b'};
	
	for i=1:numel(Tparams_fields)
		Tparams.(Tparams_fields{i}) = params.(Tparams_fields{i});
	end

	Tparams.h = h;
	Tparams.yn = yn;
	Tparams.zn = zn;
	Tparams.Ybar = Ybar;
	Tparams.Z = Z;
	Tparams.Utilde = Utilde;

	T = @(x) T1_lobattoIIIAIIIB(x,Tparams);
	DT = @(x) DT1_lobattoIIIAIIIB(x,Tparams);

	if strcmpi(iterative_scheme, 'newton')
		[root iter] = newton(T,DT,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton')
		[root iter] = simplified_newton(T,DT(x0),x0,tol);
	elseif strcmpi(iterative_scheme, 'kelley')
		combined = @(x) combined1(x,T,DT);
		tol_nsold = [1.d-10,1.d-10];
		params_nsold = [40, -1, 0.5, 0];
		%params_nsold = [40, -1, 0.5, 1];
		[root, it_hist, ierr, x_hist] = nsold(x0,combined,tol_nsold,params_nsold);
		iter = size(it_hist,1);
	else
		error('Error (solve_second_system_lobattoIIIAIIIB): Invalid iterative_scheme')
	end

	% unpack root 
	znp = root(1:nz);
	Us = root(nz+1:nz+nu);
end

