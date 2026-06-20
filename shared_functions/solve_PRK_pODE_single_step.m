function [ynp znp Ynp Znp iter] = solve_PRK_pODE_single_step(params, h, yn, zn, x0, iterative_scheme)
	%function [ynp znp Ynp Znp iter] = solve_PRK_pODE_single_step(params, h, yn, zn, x0, iterative_scheme)
	%
	% Solves the partitioned ODE y'=f(y,z), z'=g(y,z) by a PRK method
	% over a single step.
	% 
	% This is meant to work with any s-stage implicit PRK method.  	
	%
	% Inputs: 
	%	params: Structured array of parameters.
	%	h: step size.
	%	yn, zn: The most recent step.  
	%	x0: The starting approximation for the iterative
	%		solver for the system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton' or 'simplified newton 1': Uses the
	%			matrix DF(x0)
	%
	% Outputs:
	%	ynp, znp: The update step
	% 	Ynp, Znp: Vector of the internal stages
	% 	iter: Number of iterations used by the iterative scheme.

	% unpack params
	problem_str = params.problem_str;
	f = params.f;
	g = params.g;
	dfdy = params.dfdy;
	dfdz = params.dfdz;
	dgdy = params.dgdy;
	dgdz = params.dgdz;
	ny = params.ny;
	nz = params.nz;
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
	% The argument of H is [Y;Z], where Y=[Y1;...;Ys] and Z=[Z1;...;Zs].

	% subset of params except for h, yn, and zn
	Hparams_fields = {'ny', 'nz', 'A', 'Ahat', 'f', 'g', 'dfdy', 'dfdz', 'dgdy', 'dgdz'};
	for i=1:numel(Hparams_fields)
		Hparams.(Hparams_fields{i}) = params.(Hparams_fields{i});
	end

	Hparams.h = h;
	Hparams.yn = yn;
	Hparams.zn = zn;

	if strcmpi(iterative_scheme, 'newton')
		H = @(x) H1_pODE(x,Hparams);
		DH = @(x) DH1_pODE(x,Hparams);
		[root iter] = newton(H,DH,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton') || strcmpi(iterative_scheme, 'simplified newton 1')
		H = @(x) H1_pODE(x,Hparams);
		DH = @(x) DH1_pODE(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);
	elseif strcmpi(iterative_scheme, 'kelley')
		H = @(x) H1_pODE(x,Hparams);
		DH = @(x) DH1_pODE(x,Hparams);

		combined = @(x) combined1(x,H,DH);
		tol_nsold = [1.d-10,1.d-10];
		params_nsold = [40, -1, 0.5, 0];
		%params_nsold = [40, -1, 0.5, 1];
		[root, it_hist, ierr, x_hist] = nsold(x0,combined,tol_nsold,params_nsold);
		iter = size(it_hist,1);
	else
		error('Error (solve_PRK_pODE_single_step): Unknown iterative scheme');
	end

	Ynp = root(1:ny*s);
	Znp = root(ny*s+1:(ny+nz)*s);

	% Find the update step.
	% Note that computeF can be used for both partitioned ODEs and 
	% index 3 DAEs.
	ynp = yn + h*kron(b',eye(ny))*computeF(Ynp,Znp,f,ny,nz,s);
	znp = zn + h*kron(b',eye(nz))*computeF(Ynp,Znp,g,ny,nz,s);
end


