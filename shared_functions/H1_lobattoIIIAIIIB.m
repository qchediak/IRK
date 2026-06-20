function v = H1_lobattoIIIAIIIB(x, Hparams)
	% function v = H1_lobattoIIIAIIIB(x, Hparams)
	%
	% RHS of first nonlinear system 0=H1(x) of Lobatto IIIA-IIIB applied
	% to an index 3 DAE.
	%
	% Inputs:
	%	x = (Ybar, Z, Utilde): The argument of H, where
	%		Ybar = [Y2; ...; Ys], 
	%		Z = [Z1; ...; Zs], 
	% 		Utilde = [U1; ...; U_{s-1}].
	%
	%	Hparams: Structured array with fields listed below.
	%
	% Ouputs:
	%	v: The value of H1(x).

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

	% unpack x
	Ybar = x(1:(s-1)*ny);
	Z = x((s-1)*ny+1:(s-1)*ny+s*nz);
	Utilde = x((s-1)*ny+s*nz+1:(s-1)*ny+s*nz+(s-1)*nu);

	Y1 = yn; % previous step
	Y = [Y1; Ybar];

	% make a stand-in for Us for the computeK function
	Us_nan = nan(nu,1);

	F = computeF(Y,Z,f,ny,nz,s);
	K_YZU = computeK(Y,Z,[Utilde; Us_nan],k,ny,nz,nu,s);
	Ktilde = K_YZU(1:(s-1)*nz); % only keep the first s-1 sub-vectors
	G_Y = computeG(Y,g,ny,s); % this is [g(Y1);...;g(Ys)]
	G = G_Y(nu+1:s*nu); % only keep [g(Y2);...;g(Ys)]

	v1 = Ybar - kron(ones(s-1,1),yn) - h*kron(A0,eye(ny))*F;
	v2 = Z - kron(ones(s,1),zn) - h*kron(A0hat,eye(nz))*Ktilde;
	v3 = (1/h^2)*G; % Note the constratin equation has been scaled

	v = [v1; v2; v3];
end


