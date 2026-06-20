function v = H1_pODE(x,Hparams)
	% function v = H1_pODE(x,Hparams)
	%
	% The nonlinear system to solve a partitioned ODE by a PRK method.
	%
	% Note that the nonlinear system changes with each step, because yn
	% and zn are updated.  
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, h, yn, zn, A, Ahat, f, g, dfdy, dfdz, dgdy, dgdz, 
	%	where: 
	%			ny, nz: dimensions of y and z
	%			h: stepsize 
	%			yn, zn: the previous step
	% 			A, Ahat: the PRK coefficient matrices
	% 			f, g: RHS functions of the partitioned ODE 
	%				y'=f(y,z), z'=g(y,z)
	%			dfdy, dfdz, dgdy, dgdz: partial derivatives
	%
	% 	x: The argument of H.  Here x=[Y;Z], where Y=[Y1;...;Ys] and
	% 		Z=[Z1;...;Zs].
	%
	% Outputs:
	%	v: The value H(x).

	% unpack Hparams
	ny = Hparams.ny;
	nz = Hparams.nz;
	h = Hparams.h;
	yn = Hparams.yn;
	zn = Hparams.zn;
	A = Hparams.A;
	Ahat = Hparams.Ahat;
	f = Hparams.f;
	g = Hparams.g;
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dgdy = Hparams.dgdy;
	dgdz = Hparams.dgdz;
		
	s = size(A,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*(ny+nz)
		error('Error (H1_pODE): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (H1_pODE): x has too many columns.')
	end

	Y = x(1:ny*s);
	Z = x(ny*s+1:(ny+nz)*s);

	% calculate H(x) with above parameters
	v1 = Y - kron(ones(s,1), yn) - h * kron(A, eye(ny)) * computeF(Y,Z,f,ny,nz,s);
	v2 = Z - kron(ones(s,1), zn) - h * kron(Ahat, eye(nz)) * computeF(Y,Z,g,ny,nz,s);

	v = [v1; v2];
end

