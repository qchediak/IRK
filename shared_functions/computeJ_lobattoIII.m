function J = computeJ_lobattoIII(Hparams)
	% function J = computeJ_lobattoIII(Hparams)
	% Finds the matrix J for use in the root finder. 
	% Used with the projected Lobatto III method.
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, nu, h, yn, zn, un, M, Mhat, f, k, g, dfdy,
	% 		dfdz, dkdy, dkdz, dkdu, and dgdy, 
	% 	where 
	%		ny,nz,nu: dimensions of y,z,u
	% 		h: stepsize
	%		(yn,zn,un): the previous step
	%		M, Mhat: Matrices made from PRK coefficients
	%		f, k, g, derivatives: functions from the index 3 DAE
	%
	% Outputs:
	%	J: The matrix J used in the iterative root finder.
	%		An approximation to the Jacobian DH(x).

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

	% make J as a 3x3 block matrix

	J11 = eye(ny*s);
	J12 = -h*kron(M,dfdz(yn,zn));
	J13 = zeros(ny*s,nu*s);

	J21 = -h*kron(Mhat*diag(ones(s-2,1),-2),dkdy(yn,zn,un));
	J22 = eye(nz*s);
	J23 = -h*kron(Mhat,dkdu(yn,zn,un));

	J31 = kron(M*diag(ones(s-1,1),-1), dgdy(yn)*dfdy(yn,zn));
	J32 = kron(M, dgdy(yn)*dfdz(yn,zn));
	J33 = zeros(nu*s);

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];

end




