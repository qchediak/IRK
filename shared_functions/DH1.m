function J = DH1(x,Hparams)
	% function J = DH1(x,Hparams)
	% Jacobian of the function H1.
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, nu, h, yn, zn, un, A, Ahat, f, k, g, dfdy, dfdz,
	%		dkdy, dkdz, dkdu, and dgdy, where: 
	%			ny,nz,nu: dimensions of	y,z,u 
	%			h: stepsize 
	%			(yn,zn,un): the previous step,
	% 			A, Ahat: the PRK coefficient matrices
	% 			f, k, g and derivatives: Functions the index 3 problem.
	%
	% 	x: The argument of DH. Here x = [Y; Z; U].
	%
	% Outputs:
	%	J: The Jacobian matrix of the function H1.

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
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dkdy = Hparams.dkdy;
	dkdz = Hparams.dkdz;
	dkdu = Hparams.dkdu;
	dgdy = Hparams.dgdy;

	s = size(A,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*(ny+nz+nu)
		error('Error (DH1): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DH1): x has too many columns.')
	end

	% unpack x
	Y = x(1:ny*s);				% for s=2, ny=nz=2, nu=1, this is 1:4
	Z = x(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
	U = x(ny*s+nz*s+1:end);		% for s=2, ny=nz=2, nu=1, this is 9:10

	% compute the derivatives 
	dFdY = compute_dFdY(Y,Z,dfdy,ny,nz,s);
	dFdZ = compute_dFdZ(Y,Z,dfdz,ny,nz,s);

	dKdY = compute_dKdY(Y,Z,U,dkdy,ny,nz,nu,s);
	dKdZ = compute_dKdZ(Y,Z,U,dkdz,ny,nz,nu,s);
	dKdU = compute_dKdU(Y,Z,U,dkdu,ny,nz,nu,s);

	dGdY = compute_dGdY(Y,dgdy,ny,s);

	% calculate the Jacobian as a 3x3 block matrix
	J11 = eye(ny*s) - h * kron(A,eye(ny)) * dFdY;
	J12 = -h * kron(A, eye(ny)) * dFdZ;
	J13 = zeros(ny*s,nu*s);

	J21 = -h * kron(Ahat, eye(nz)) * dKdY;
	J22 = eye(nz*s) - h * kron(Ahat,eye(nz)) * dKdZ;
	J23 = -h * kron(Ahat, eye(nz)) * dKdU;

	J31 = dGdY;
	J32 = zeros(nu*s,nz*s);
	J33 = zeros(nu*s,nu*s);

	J = [J11 J12 J13; 
		J21 J22 J23; 
		J31 J32 J33];
end


