function J = DT1_lobattoIIIAIIIB(x,Tparams)
	% function J = DT1_lobattoIIIAIIIB(x,Tparams)
	% Jacobian of T1_lobattoIIIAIIIB.
	%
	% Inputs:
	%	x = (znp, Us): The arguemnt of DT.
	%	Tparams: Structured array with the following fields:
	%		ny, nz, nu: dimensions of y, z, u
	%		h: stepsize
	%		yn, zn: The previous step
	%		f, k, g: Functions of the DAE
	%		dfdy, dfdz, dkdy, dkdz, dkdu, dgdy: Derivatives of f,k,g
	%		Ybar, Z, Utilde: Solution to previous nonlinear system.
	%
	%	Outputs:
	%		J: The value DT(x).

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

	Ys = Y((s-1)*ny+1:end);
	Zs = Z((s-1)*nz+1:end);
	ynp = Ys;

	% compute the Jacobian as a 2x2 block matrix
	J11 = eye(nz);
	J12 = -h*b(s)*dkdu(Ys,Zs,Us);
	J21 = (1/h)*dgdy(ynp)*dfdz(ynp,znp);
	J22 = zeros(nu,nu);

	J = [J11 J12; J21 J22];
end

