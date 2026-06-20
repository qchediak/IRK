function v = H1_ODE(x,Hparams)
	% function v = H1_ODE(x,Hparams)
	%
	% The nonlinear system to solve an ODE by an IRK method.
	%
	% Note that the nonlinear system changes with each step, because yn
	% is updated.  
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, h, yn, A, f, dfdy, where: 
	%			ny: dimension of y
	%			h: stepsize 
	%			yn: the previous step,
	% 			A: the IRK coefficient matrix
	% 			f: RHS function of the ODE y'=f(y)
	%			dfdy: Jacobian of f
	%
	% 	x: The argument of H. Here x=Y=[Y1;....;Ys].
	%
	% Outputs:
	%	v: The value H(x).

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
		error('Error (H1_ODE): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (H1_ODE): x has too many columns.')
	end

	Y = x;

	% calculate H(x) with above parameters
	v = Y - kron(ones(s,1), yn) - h * kron(A, eye(ny)) * computeF_ODE(Y,f,ny,s);
end

