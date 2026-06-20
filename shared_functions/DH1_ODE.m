function J = DH1_ODE(x,Hparams)
	% function J = DH1_ODE(x,Hparams)
	% Jacobian of the function H1.
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, h, yn, A, f, dfdy, where: 
	%			ny: dimension of y
	%			h: stepsize 
	%			yn: the previous step,
	% 			A: the IRK coefficient matrix
	% 			f: the RHS function of y'=f(y)
	%			dfdy: The Jacobian of f.
	%
	% 	x: The argument of DH. Here x=Y=[Y1;...;Ys].
	%
	% Outputs:
	%	J: The Jacobian matrix of the function H1_ODE.

	% unpack Hparams
	ny = Hparams.ny;
	h = Hparams.h;
	yn = Hparams.yn;
	A = Hparams.A;
	f = Hparams.f;
	dfdy = Hparams.dfdy;

	s = size(A,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*ny
		error('Error (DH1_ODE): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DH1_ODE): x has too many columns.')
	end

	Y = x;	

	% compute the derivatives 
	dFdY = compute_dFdY_ODE(Y,dfdy,ny,s);

	% calculate the Jacobian 
	J = eye(ny*s) - h * kron(A,eye(ny)) * dFdY;
end


