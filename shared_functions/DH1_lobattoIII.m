function J = DH1_lobattoIII(x, Hparams)
	% function J = DH1_lobattoIII(x, Hparams)
	%
	% Jacobian of the function H1_lobattoIII. The function H1_lobattoIII
	% has no scaling of the constraint, so this Jacobian does not
	% either.
	%
	% Note that the nonlinear system changes with each step, because
	% (yn,zn) are updated.  This is why we need to pass Hparams.
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
	%	J: The value DH(x).

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
		error('Error (DH1_lobattoIII): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DH1_lobattoIII): x has too many columns.')
	end

	% check that s \geq 3
	if s <= 2
		error('Error (DH1_lobattoIII): This method requires s >= 3')
	end

	% unpack x
	Yhat = x(1:ny*s);		
	Z = x(ny*s+1:ny*s+nz*s);	
	Uhat = x(ny*s+nz*s+1:ny*s+nz*s+nu*s);

	% define the "zero stage" and full vector Y
	Y0 = yn; % Y0 = Y1 = yn
	Z0 = Z(1:nz); % Z0 = Z1
	U0 = Uhat(1:nu); 
	Y = [yn; Yhat(1:(s-1)*ny)]; % [Y1; ... Ys]

	% --------------------
	% make the derivatives
	% --------------------

	% dFhat/dYhat
	Y_part_dFdY = Yhat(1:(s-1)*ny); % [Y2; ...; Ys]
	Z_part_dFdY = Z(nz+1:s*nz); % [Z2; ...; Zs]
	
	dFhatdYhat_part = compute_dFdY(Y_part_dFdY, Z_part_dFdY, dfdy, ny, nz, s-1);
	dFhatdYhat = [zeros(ny, (s-1)*ny) zeros(ny);
					dFhatdYhat_part zeros((s-1)*ny, ny)];

	% dFhat/dZ
	dFhatdZ = compute_dFdZ(Y,Z,dfdz,ny,nz,s);

	% dKhat/dYhat
	Y_part_dKdY = Y(ny+1:(s-1)*ny); % [Y_2; ...; Y_{s-1}]
	Z_part_dKdY = Z(nz+1:(s-1)*nz); % [Z_2; ...; Z_{s-1}]
	U_part_dKdY = Uhat(2*nu+1:s*nu); % [U_2; ...; U_{s-1}]

	dKhatdYhat_part = compute_dKdY(Y_part_dKdY, Z_part_dKdY, U_part_dKdY, dkdy, ny, nz, nu, s-2);

	dKhatdYhat = [zeros(2*nz,(s-2)*ny) zeros(2*nz,2*ny);
					dKhatdYhat_part zeros((s-2)*m,2*ny)];

	% dKhat/dZ
	dKhatdZ_11 = [dkdz(Y0,Z0,U0) zeros(nz, (s-2)*nz)];

	Y_part_dKdZ = Y(1:(s-1)*ny); % [Y_1; ...; Y_{s-1}]
	Z_part_dKdZ = Z(1:(s-1)*nz); % [Z_1; ...; Z_{s-1}]
	U_part_dKdZ = Uhat(nu+1:nu*s); % [U_1; ...; U_{s-1}]

	dKhatdZ_21 = compute_dKdZ(Y_part_dKdZ, Z_part_dKdZ, U_part_dKdZ, dkdz, ny, nz, nu, s-1);

	dKhatdZ = [dKhatdZ_11 zeros(nz); dKhatdZ_21 zeros((s-1)*nz,nz)];

	% dKhat/dUhat
	Y_part_dKdU = Y(1:(s-1)*ny); % [Y_1; ...; Y_{s-1}]
	Y_part_dKdU = [Y0; Y_part_dKdU]; % [Y_0; ...; Y_{s-1}]

	Z_part_dKdU = Z(1:(s-1)*nz); % [Z_1; ...; Z_{s-1}]
	Z_part_dKdU = [Z0; Z_part_dKdU]; % [Z_0; ...; Z_{s-1}]

	dKhatdUhat = compute_dKdU(Y_part_dKdU, Z_part_dKdU, Uhat, ny, nz, nu, s);

	% dG/dYhat
	dGdYhat = computeG(Yhat, g, ny, s);

	% --------------------------------------------
	% calculate the Jacobian as a 3x3 block matrix
	% --------------------------------------------
	J11 = eye(ny*s) - h * kron(M,eye(ny)) * dFhatdYhat;
	J12 = -h * kron(M, eye(ny)) * dFhatdZ;
	J13 = zeros(ny*s,nu*s);

	J21 = -h * kron(Mhat, eye(nz)) * dKhatdYhat;
	J22 = eye(nz*s) - h * kron(Mhat,eye(nz)) * dKhatdZ;
	J23 = -h * kron(Mhat, eye(nz)) * dKhatdUhat;

	J31 = dGdYhat;
	J32 = zeros(nu*s,nz*s);
	J33 = zeros(nu*s);

	J = [J11 J12 J13; 
		J21 J22 J23; 
		J31 J32 J33];
end


