function v = T1_lobattoIIIAIIIB(x,Tparams)
	% function v = T1_lobattoIIIAIIIB(x,Tparams)
	% RHS of second nonlinear system 0=T(x) to solve for Lobatto IIIA-IIIB.
	%
	% Inputs:
	%	x = (znp, Us): The arguemnt of T.
	%	Tparams: Structured array with the following fields:
	%		ny, nz, nu: dimensions of y, z, u
	%		h: stepsize
	%		yn, zn: The previous step
	%		f, k, g: Functions of the DAE
	%		dfdy, dfdz, dkdy, dkdz, dkdu, dgdy: Derivatives of f,k,g
	%		Ybar, Z, Utilde: Solution to previous nonlinear system.
	%		b: Vector of weights.
	%
	%	Outputs:
	%		v: The value T(x).

	% unpack Tparams
	ny = Tparams.ny;
	nz = Tparams.nz;
	nu = Tparams.nu;
	h = Tparams.h;
	yn = Tparams.yn;
	zn = Tparams.zn;
	f = Tparams.f;
	k = Tparams.k;
	g = Tparams.g;
	dfdy = Tparams.dfdy;
	dfdz = Tparams.dfdz;
	dkdy = Tparams.dkdy;
	dkdz = Tparams.dkdz;
	dkdu = Tparams.dkdu;
	dgdy = Tparams.dgdy;
	Ybar = Tparams.Ybar;
	Z = Tparams.Z;
	Utilde = Tparams.Utilde;
	b = Tparams.b;

	s = length(b);

	% unpack x
	znp = x(1:nz);
	Us = x(nz+1:nz+nu);

	Y1 = yn; % previous step
	Y = [Y1; Ybar];
	U = [Utilde; Us];

	ynp = Y((s-1)*ny+1:end); % ynp = Ys

	% Find T(x)
	v1 = znp - zn - h*kron(b',eye(nz))*computeK(Y,Z,U,k,ny,nz,nu,s);
	v2 = (1/h)*dgdy(ynp)*f(ynp,znp);

	v = [v1; v2];
end

