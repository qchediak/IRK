function J = DH1_pODE(x,Hparams)
	%function J = DH1_pODE(x,Hparams)
	% Jacobian of the function H1_pODE.
	%
	% Inputs:
	%	Hparams: Structured array with fields
	%		ny, nz, h, yn, zn, A, Ahat, f, g, dfdy, dfdz, dgdy, dgdz, 
	%	where: 
	%			ny, nz: dimensions of y and z
	%			h: stepsize 
	%			yn, zn: the previous step
	% 			A, Ahat: the PRK coefficient matrices
	% 			f, g: RHS functions of the partitioned ODE 
	%				y'=f(y,z), z'=g(y,z)
	%			dfdy, dfdz, dgdy, dgdz: partial derivatives
	%
	% 	x: The argument of H.  Here x=[Y;Z], where Y=[Y1;...;Ys] and
	% 		Z=[Z1;...;Zs].
	%
	% Outputs:
	%	J: The Jacobian matrix.

	% unpack Hparams
	ny = Hparams.ny;
	nz = Hparams.nz;
	h = Hparams.h;
	yn = Hparams.yn;
	zn = Hparams.zn;
	A = Hparams.A;
	Ahat = Hparams.Ahat;
	f = Hparams.f;
	g = Hparams.g;
	dfdy = Hparams.dfdy;
	dfdz = Hparams.dfdz;
	dgdy = Hparams.dgdy;
	dgdz = Hparams.dgdz;

	s = size(A,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= s*(ny+nz)
		error('Error (DH1_pODE): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DH1_ODE): x has too many columns.')
	end

	Y = x(1:s*ny);	
	Z = x(s*ny+1:s*(ny+nz));	

	% Compute the derivatives. 
	% Note that compute_dFdY can be used for both partitioned ODEs and
	% index 3 DAEs.
	dFdY = compute_dFdY(Y,Z,dfdy,ny,nz,s);
	dFdZ = compute_dFdY(Y,Z,dfdz,ny,nz,s);
	dGdY = compute_dFdY(Y,Z,dgdy,ny,nz,s);
	dGdZ = compute_dFdY(Y,Z,dgdz,ny,nz,s);

	% calculate the Jacobian as a 2x2 block matrix
	J11 = eye(ny*s) - h * kron(A,eye(ny)) * dFdY;
	J12 = -h * kron(A,eye(ny)) * dFdZ;

	J21 = -h * kron(Ahat,eye(nz)) * dGdY;
	J22 = eye(nz*s) - h * kron(Ahat,eye(nz)) * dGdZ;

	J = [J11 J12; J21 J22];
end


