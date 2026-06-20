function [Ybar Z Utilde iter] = solve_first_system_lobattoIIIAIIIB(params, h, yn, zn, un, x0, iterative_scheme)
	% function [Ybar Z Utilde] = solve_first_system_lobattoIIIAIIIB(params, h, yn, zn, un, x0, iterative_scheme)
	%
	% Solve the first nonlinear system of Lobatto IIIA-IIIB.
	%
	% Inputs:
	%	params: Structured array with the following fields:
	%		h: step size.
	%		yn, zn, un: The most recent step. 	
	%		x0 = (Ybar0,Z0,Utilde0): The starting approximation for the
	%			iterative solver for the system of nonlinear equations.
	%		iterative_scheme: The iterative scheme used to solve the
	%			nonlinear system.  Can take the following values: 
	%				'newton': The standard Newton method.
	%				'simplified newton': Uses the matrix DF(x0)
	%				'simplified newton 2': Uses approximate Jacobian J
	%
	% Outputs:
	%	Ybar = (Y2,...Ys): Internal stages of the y-component
	%	Z = (Z1,...,Zs): Internal stages of the z-component
	%	Utilde = (U1,...,U_{s-1}): Internal stages of the u-component
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

	A0 = A(2:s,:);
	A0hat = Ahat(:,1:s-1);

	% subset of params except for A0, A0hat, h, yn, zn, un
	Hparams_fields = {'ny', 'nz', 'nu', 'f', 'k', 'g', 'dfdy', 'dfdz', 'dkdy', 'dkdz', 'dkdu', 'dgdy'};
	for i=1:numel(Hparams_fields)
		Hparams.(Hparams_fields{i}) = params.(Hparams_fields{i});
	end

	Hparams.A0 = A0;
	Hparams.A0hat = A0hat;
	Hparams.h = h;
	Hparams.yn = yn;
	Hparams.zn = zn;
	Hparams.un = un;

	if strcmpi(iterative_scheme, 'newton')

		H = @(x) H1_lobattoIIIAIIIB(x,Hparams);
		DH = @(x) DH1_lobattoIIIAIIIB(x,Hparams);
		[root iter] = newton(H,DH,x0,tol);

	elseif strcmpi(iterative_scheme, 'simplified newton 0') 

		H = @(x) H0_lobattoIIIAIIIB(x,Hparams);
		DH = @(x) DH0_lobattoIIIAIIIB(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);

	elseif strcmpi(iterative_scheme, 'simplified newton 1') || strcmpi(iterative_scheme, 'simplified newton')

		H = @(x) H1_lobattoIIIAIIIB(x,Hparams);
		DH = @(x) DH1_lobattoIIIAIIIB(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);

	elseif strcmpi(iterative_scheme, 'simplified newton 2')

		H = @(x) H2_lobattoIIIAIIIB(x,Hparams);
		J = computeJ_lobattoIIIAIIIB(Hparams);
		[root iter] = simplified_newton(H,J,x0,tol);

	elseif strcmpi(iterative_scheme, 'kelley')

		H = @(x) H1_lobattoIIIAIIIB(x,Hparams);
		DH = @(x) DH1_lobattoIIIAIIIB(x,Hparams);

		combined = @(x) combined1(x,H,DH);
		tol_nsold = [1.d-10,1.d-10];
		params_nsold = [40, -1, 0.5, 0];
		%params_nsold = [40, -1, 0.5, 1];
		[root, it_hist, ierr, x_hist] = nsold(x0,combined,tol_nsold,params_nsold);
		iter = size(it_hist,1);

	else

		error('Error (solve_first_system_lobattoIIIAIIIB): Invalid iterative_scheme.')

	end

	% unpack root
	Ybar = root(1:(s-1)*ny);
	Z = root((s-1)*ny+1:(s-1)*ny+s*nz);
	Utilde = root((s-1)*ny+s*nz+1:(s-1)*ny+s*nz+(s-1)*nu);
end
	

