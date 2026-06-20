function [ynp Ynp iter] = solve_IRK_ODE_single_step(params, h, yn, x0, iterative_scheme)
	%function [ynp Ynp iter] = solve_IRK_ODE_single_step(params, h, yn, x0, iterative_scheme)
	%
	% Solves the ODE y'=f(y) by an IRK method over a single step.
	% 
	% This is meant to work with any s-stage implicit RK method.  	
	%
	% Inputs: 
	%	params: Structured array of parameters.
	%	h: step size.
	%	yn: The most recent step.  
	%	x0: The starting approximation for the iterative
	%		solver for the system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton' or 'simplified newton 1': Uses the
	%			matrix DF(x0)
	%
	% Outputs:
	%	ynp: The update step
	% 	Ynp: Vector of the internal stages
	% 	iter: Number of iterations used by the iterative scheme.

	% unpack params
	problem_str = params.problem_str;
	f = params.f;
	dfdy = params.dfdy;
	ny = params.ny;
	method_str = params.method_str;
	b = params.b;
	c = params.c;
	A = params.A;
	tol = params.tol;

	s = size(A,1);

	% ---------------------------
	% Solve the nonlinear system
	% ---------------------------

	% The function H needs to be redefined at every step using Hparams.
	% The argument of H is Y, where Y = [Y1;...;Ys].

	% subset of params except for h and yn
	Hparams_fields = {'ny', 'A', 'f', 'dfdy'};
	for i=1:numel(Hparams_fields)
		Hparams.(Hparams_fields{i}) = params.(Hparams_fields{i});
	end

	Hparams.h = h;
	Hparams.yn = yn;

	if strcmpi(iterative_scheme, 'newton')
		H = @(x) H1_ODE(x,Hparams);
		DH = @(x) DH1_ODE(x,Hparams);
		[root iter] = newton(H,DH,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton') || strcmpi(iterative_scheme, 'simplified newton 1')
		H = @(x) H1_ODE(x,Hparams);
		DH = @(x) DH1_ODE(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);
	else
		error('Error (solve_IRK_ODE_single_step): Unknown iterative scheme');
	end

	Ynp = root;

	% Find the update step
	ynp = yn + h*kron(b',eye(ny))*computeF_ODE(Ynp,f,ny,s);
end


