function J = computeJ3(Hparams)
	% function J = computeJ3(Hparams)
	% Finds the matrix J for use in the root finder. 
	% This is like computeJ, except the returned matrix has some entries
	% set to zero as a further approximation.
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, nu, h, yn, zn, un, A, Ahat, f, k, g, dfdy,
	% 		dfdz, dkdy, dkdz, dkdu, and dgdy, 
	% 	where 
	%		ny,nz,nu: dimensions of y,z,u
	% 		h: stepsize
	%		(yn,zn,un): the previous step
	%		A, Ahat: PRK coefficient matrices
	%		f, k, g, derivatives: functions from the index 3 DAE
	%
	% Outputs:
	%	J: The matrix J used in the iterative root finder.

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
	dfdz = Hparams.dfdz;
	dkdz = Hparams.dkdz;
	dkdu = Hparams.dkdu;
	dgdy = Hparams.dgdy;

	s = size(A,1);

	% compute J as a 3x3 block matrix

	J11 = eye(ny*s);
	J12 = -h*kron(A,dfdz(yn,zn));
	J13 = zeros(ny*s,nu*s);

	J21 = zeros(s*nz,s*ny);
	J22 = eye(nz*s);
	J23 = -h*kron(Ahat,dkdu(yn,zn,un));

	J31 = zeros(s*nu,s*ny);
	J32 = kron(A, dgdy(yn)*dfdz(yn,zn));
	J33 = zeros(nu*s);

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];

end




