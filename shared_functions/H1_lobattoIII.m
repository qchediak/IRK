function v = H1(x,Hparams)
	% function v = H1(x,Hparams)
	%
	% The RHS of the nonlinear system to solve for the projected 
	% Lobatto III method applied to an index 3 DAE. In this nonlinear
	% system, the constraint equation is NOT scaled.
	%
	% Note that the nonlinear system changes with each step, because
	% (yn,zn) are updated.  
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, nu, h, yn, zn, un, M, Mhat, f, k, g, dfdy, dfdz,
	%		dkdy, dkdz, dkdu, and dgdy, where: 
	%			ny,nz,nu: dimensions of	y,z,u 
	%			h: stepsize 
	%			(yn,zn,un): the previous step,
	% 			M, Mhat: Matrices made from PRK coefficients
	% 			f, k, g and derivatives: Functions the index 3 problem.
	%
	% 	x: The argument of H. Here x = [Yhat; Z; Uhat], where
	%		Yhat=(Y2,...,Ys,ynp), Z=(Z1,...,Zs), and
	%		Uhat = (U0, U1,...,U_{s-1}).
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
	M = Hparams.M;
	Mhat = Hparams.Mhat;
	f = Hparams.f;
	k = Hparams.k;
	g = Hparams.g;
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dkdy = Hparams.dkdy;
	dkdz = Hparams.dkdz;
	dkdu = Hparams.dkdu;
	dgdy = Hparams.dgdy;

	s = size(M,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*(ny+nz+nu)
		error('Error (H1_lobattoIII): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (H1_lobattoIII): x has too many columns.')
	end

	% unpack x
	Yhat = x(1:ny*s);		
	Z = x(ny*s+1:ny*s+nz*s);	
	Uhat = x(ny*s+nz*s+1:ny*s+nz*s+nu*s);

	Fhat = computeFhat(Yhat,Z,f,ny,nz,s);
	Khat = computeKhat(Yhat,Z,Uhat,k,ny,nz,nu,s);
	G = computeG(Yhat,g,ny,s);

	% calculate H(x) with above parameters
	v1 = Yhat - kron(ones(s,1), yn) - h * kron(M, eye(ny)) * Fhat;
	v2 = Z - kron(ones(s,1), zn) - h * kron(Mhat, eye(nz)) * Khat;
	v3 = G;

	v = [v1; v2; v3];
end

