function [ynp znp Ynp Znp Unp iter] = solve_PRK_DAE_single_step(params, h, yn, zn, un, x0, iterative_scheme)
	% function [ynp znp Ynp Znp Unp] = solve_IRK_DAE_single_step(params, h, yn, zn, x0, iterative_scheme)
	%
	% Solves the index 3 DAE y'=f(y,z), z'=k(y,z,u), 0=g(y).
	% 
	% This is meant to work with any s-stage PRK method.  However, it
	% does NOT exploit the special structure of the Lobatto IIIA-IIIB
	% method, so using this function with Lobatto IIIA-IIIB is
	% suboptimal.
	%
	% Inputs: 
	%	params: Structured array of parameters.
	%	h: step size.
	%	(yn, zn, un): The most recent step.  
	%	x0 = (Y0,Z0,U0): The starting approximation for the iterative
	%		solver for the system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton' or 'simplified newton 1': Uses the
	%			matrix DF(x0)
	%		'simplified newton 2': Uses a matrix J which exploits the
	%			structure of the index 3 problem and approximates h as small.
	%
	% Outputs:
	%	ynp, znp: The update step
	% 	Ynp, Znp, Unp: Vectors of the internal stages
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

	% The function H needs to be redefined at every step using Hparams.
	% The argument of H is x = [Y; Z; U],
	% where Y = [Y1;...;Ys], Z = [Z1;...;Zs], U = [U1;...;Us].

	% subset of params except for h, yn, zn, un
	Hparams_fields = {'ny', 'nz', 'nu', 'A', 'Ahat', 'f', 'k', 'g', 'dfdy', 'dfdz', 'dkdy', 'dkdz', 'dkdu', 'dgdy'};
	for i=1:numel(Hparams_fields)
		Hparams.(Hparams_fields{i}) = params.(Hparams_fields{i});
	end

	Hparams.h = h;
	Hparams.yn = yn;
	Hparams.zn = zn;
	Hparams.un = un;

	if strcmpi(iterative_scheme, 'newton')
		H = @(x) H1(x,Hparams);
		DH = @(x) DH1(x,Hparams);
		[root iter] = newton(H,DH,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton') || strcmpi(iterative_scheme, 'simplified newton 1')
		H = @(x) H1(x,Hparams);
		DH = @(x) DH1(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton 2') 
		H = @(x) H2(x,Hparams); % alternate formulation of the nonlinear system
		J = computeJ(Hparams); % can be used in simplified Newton
		[root iter] = simplified_newton(H,J,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton 3') 
		H = @(x) H2(x,Hparams); % alternate formulation of the nonlinear system
		J = computeJ3(Hparams); 
		[root iter] = simplified_newton(H,J,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton index 3') 
		H = @(x) H2(x,Hparams); % alternate formulation of the nonlinear system
		[root iter] = simplified_newton_index3(H,x0,tol,Hparams);
	elseif strcmpi(iterative_scheme, 'kelley')
		H = @(x) H1(x,Hparams);
		DH = @(x) DH1(x,Hparams);

		combined = @(x) combined1(x,H,DH);
		tol_nsold = [1.d-10,1.d-10];
		params_nsold = [40, -1, 0.5, 0];
		%params_nsold = [40, -1, 0.5, 1];
		[root, it_hist, ierr, x_hist] = nsold(x0,combined,tol_nsold,params_nsold);
		iter = size(it_hist,1);
	else
		error('Error (solve_PRK_DAE_single_step): Unknown iterative scheme');
	end

	% unpack root
	Ynp = root(1:ny*s);			% for s=2, ny=nz=2, nu=1, this is 1:4
	Znp = root(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
	Unp = root(ny*s+nz*s+1:end);	% for s=2, ny=nz=2, nu=1, this is 9:10

	% Find the update step
	ynp = yn + h*kron(b',eye(ny))*computeF(Ynp,Znp,f,ny,nz,s);
	znp = zn + h*kron(b',eye(nz))*computeK(Ynp,Znp,Unp,k,ny,nz,nu,s);
end

