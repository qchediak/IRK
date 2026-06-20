function v = H1(x,Hparams)
	% function v = H1(x,Hparams)
	%
	% The nonlinear system to solve for an index 3 DAE.
	% In this nonlinear system, the constraint equation
	% is NOT scaled by 1/h.
	%
	% Note that the nonlinear system changes with each step, because (yn,zn) are
	% updated.  
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, nu, h, yn, zn, un, A, Ahat, f, k, g, dfdy, dfdz,
	%		dkdy, dkdz, dkdu, and dgdy, where: 
	%			ny,nz,nu: dimensions of	y,z,u 
	%			h: stepsize 
	%			(yn,zn,un): the previous step,
	% 			A, Ahat: the PRK coefficient matrices
	% 			f, k, g and derivatives: Functions the index 3 problem.
	%
	% 	x: The argument of H. Here x = [Y; Z; U].
	%
	% Outputs:
	%	v: The value H(x).

	% unpack Hparams
	ny = Hparams.ny;
	nz = Hparams.nz;
	nu = Hparams.nu;
	h = Hparams.h;
	yn = Hparams.yn;
	zn = Hparams.zn;
	un = Hparams.un;
	A = Hparams.A;
	Ahat = Hparams.Ahat;
	f = Hparams.f;
	k = Hparams.k;
	g = Hparams.g;
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dkdy = Hparams.dkdy;
	dkdz = Hparams.dkdz;
	dkdu = Hparams.dkdu;
	dgdy = Hparams.dgdy;

	s = size(A,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*(ny+nz+nu)
		error('Error (H1): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (H1): x has too many columns.')
	end

	% unpack x
	Y = x(1:ny*s);				% for s=2, ny=nz=2, nu=1, this is 1:4
	Z = x(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
	U = x(ny*s+nz*s+1:ny*s+nz*s+nu*s);	% for s=2, ny=nz=2, nu=1, this is 9:10

	% calculate H(x) with above parameters
	v1 = Y - kron(ones(s,1), yn) - h * kron(A, eye(ny)) * computeF(Y,Z,f,ny,nz,s);
	v2 = Z - kron(ones(s,1), zn) - h * kron(Ahat, eye(nz)) * computeK(Y,Z,U,k,ny,nz,nu,s);
	v3 = computeG(Y,g,ny,s);

	v = [v1; v2; v3];
end

