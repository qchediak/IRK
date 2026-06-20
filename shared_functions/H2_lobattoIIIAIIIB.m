function v = H2_lobattoIIIAIIIB(x, Hparams)
	% function v = H2_lobattoIIIAIIIB(x, Hparams)
	%
	% RHS of first nonlinear system 0=H2(x) of Lobatto IIIA-IIIB applied
	% to an index 3 DAE.  Uses a modified constraint scaled by 1/h, 
	% which is equivalent to the original constraint.
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
	%	v: The value of H(x).

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
	K_YZU = computeK(Y, Z, [Utilde; Us_nan], k, ny, nz, nu, s);
	Ktilde = K_YZU(1:(s-1)*nz); % only keep the first s-1 sub-vectors

	% first argument for computeG
	G_arg = [nan(ny,1); kron(ones(s-1,1),yn) + h * kron(A0,eye(ny))*computeF(Y,Z,f,ny,nz,s)];

	G_Y = computeG(G_arg,g,ny,s); % this is [g_1(G_arg);...;g_s(G_arg)]
	G = G_Y(nu+1:s*nu); % only keep [g_2(G_arg);...;g_s(G_arg)]

	v1 = Ybar - kron(ones(s-1,1),yn) - h*kron(A0,eye(ny))*F;
	v2 = Z - kron(ones(s,1),zn) - h*kron(A0hat,eye(nz))*Ktilde;
	v3 = (1/h)*G; % Note the constraint equation has been scaled

	v = [v1; v2; v3];
end



