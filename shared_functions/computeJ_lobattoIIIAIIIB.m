function J = computeJ_lobattoIIIAIIIB(Hparams)
	% function J = computeJ_lobattoIIIAIIIB(Hparams)
	% Finds the matrix J for use in the root finder. 
	%
	% Inputs:
	%	Hparams: Structured array with fields ny, nz, nu, h, yn, zn, un, A,
	%	Ahat, f, k, g, dfdy, dfdz, dkdy, dkdz, dkdu, and dgdy, where 
	%		ny,nz,nu: dimensions of y,z,u
	% 		h: stepsize
	%		(yn,zn,un): the previous step
	%		A, Ahat: PRK coefficient matrices
	%		f, k, g, derivatives: functions from the index 3 DAE
	%
	% Outputs:
	%	J: The approximate Jacobian J used in the iterative root finder.

	% unpack Hparams
	ny = Hparams.ny;
	nz = Hparams.nz;
	nu = Hparams.nu;
	h = Hparams.h;
	yn = Hparams.yn;
	zn = Hparams.zn;
	un = Hparams.un;
	A0 = Hparams.A0;
	A0hat = Hparams.A0hat;
	f = Hparams.f;
	k = Hparams.k;
	g = Hparams.g;
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dkdy = Hparams.dkdy;
	dkdz = Hparams.dkdz;
	dkdu = Hparams.dkdu;
	dgdy = Hparams.dgdy;

	s = size(A0hat,1);

	% make J as a 3x3 block matrix

	J11 = eye(ny*(s-1));
	J12 = -h*kron(A0,dfdz(yn,zn));
	J13 = zeros(ny*(s-1),nu*(s-1));

	J21 = -h*kron(A0hat*diag(ones(s-2,1),-1),dkdy(yn,zn,un));
	J22 = eye(nz*s);
	J23 = -h*kron(A0hat,dkdu(yn,zn,un));

	J31 = kron(A0*[zeros(1,s-1); eye(s-1)], dgdy(yn)*dfdy(yn,zn));
	J32 = kron(A0, dgdy(yn)*dfdz(yn,zn));
	J33 = zeros(nu*(s-1));

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];

end

