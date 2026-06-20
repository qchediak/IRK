function J = DH1_lobattoIIIAIIIB(x, Hparams)
	% function J = DH1_lobattoIIIAIIIB(x, Hparams)
	%
	% Jacobian of the function H1_lobattoIIIAIIIB.
	%
	% Note that the nonlinear system changes with each step, because
	% (yn,zn) are updated.  This is why we need to pass Hparams.
	%
	% Inputs:
	%	x = (Ybar, Z, Utilde): The argument of DH, where
	%		Ybar = [Y2; ...; Ys], 
	%		Z = [Z1; ...; Zs], 
	% 		Utilde = [U1; ...; U_{s-1}].
	%
	%	Hparams: Structured array with fields listed below.
	%
	% Ouputs:
	%	J: The Jacobian DH(x) of H(x).

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

	% make a stand-in for Us for the compute_dKdY etc functions
	Us_nan = nan(nu,1);
	U_standin = [Utilde; Us_nan]; 

	% compute the derivatives

	dFdY_YZ = compute_dFdY(Y,Z,dfdy,ny,nz,s);
	% keep only the lower right (s-1)x(s-1) block matrix,
	% where each block is ny x ny
	dFdY_part = dFdY_YZ(ny+1:end,ny+1:end);
	dFdY = [zeros(ny,(s-1)*ny); dFdY_part];

	dFdZ = compute_dFdZ(Y,Z,dfdz,ny,nz,s);

	dGdY_Y = compute_dGdY(Y,dgdy,ny,s);
	% keep only the lower right (s-1)x(s-1) block matrix,
	% where each block is nu x ny
	dGdY = dGdY_Y(nu+1:end,ny+1:end);

	dKdZ_YZU = compute_dKdZ(Y,Z,U_standin,dkdz,ny,nz,nu,s);
	% keep only the upper left (s-1)x(s-1) block matrix,
	% where each block is nz x nz
	dKdZ_part = dKdZ_YZU(1:(s-1)*nz,1:(s-1)*nz);

	dKdZ = [dKdZ_part zeros((s-1)*nz,nz)];

	dKdU_YZU = compute_dKdU(Y,Z,U_standin,dkdu,ny,nz,nu,s);
	% keep only the upper left (s-1)x(s-1) block matrix,
	% where each block is nz x nu
	dKdU = dKdU_YZU(1:(s-1)*nz,1:(s-1)*nu);

	if s==2
		dKdY = zeros(nz,ny);
	elseif s >= 3
		dKdY_YZU = compute_dKdY(Y,Z,U_standin,dkdy,ny,nz,nu,s);
		% keep only the middle (s-2)x(s-2) block matrix,
		% where each block is nz x ny.
		dKdY_part = dKdY_YZU(nz+1:(s-1)*nz,ny+1:(s-1)*ny);

		dKdY = [zeros(nz,(s-2)*ny) zeros(nz,ny); 
				dKdY_part zeros((s-2)*nz,ny)];
	else
		error('Error (DH1_lobattoIIIAIIIB): Invalid value of s')
	end

	% Find the Jacobian as a 3 x 3 block matrix

	J11 = eye((s-1)*ny) - h*kron(A0,eye(ny))*dFdY;
	J12 = -h*kron(A0,eye(ny))*dFdZ;
	J13 = zeros((s-1)*ny,(s-1)*nu);

	J21 = -h*kron(A0hat,eye(nz))*dKdY;
	J22 = eye(s*nz) - h*kron(A0hat,eye(nz))*dKdZ;
	J23 = -h*kron(A0hat,eye(nz))*dKdU;

	J31 = (1/h^2)*dGdY; % constraint is scaled
	J32 = zeros((s-1)*nu,s*nz);
	J33 = zeros((s-1)*nu,(s-1)*nu);

	J = [J11 J12 J13; 
		J21 J22 J23; 
		J31 J32 J33];

end

