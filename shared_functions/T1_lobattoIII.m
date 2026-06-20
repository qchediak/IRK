function v = T1_lobattoIII(x,Tparams)
	% function v = T1(x,Tparams)
	% RHS of second linear system to solve for projected Lobatto III
	% applied to an index 3 DAE.
	%
	% Note that the linear system changes with each step, because (yn,zn)
	% are updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unps]: The argument of T.
	% 	Tparams: Structured array with fields ny, nz, nu, h, yn, zn,
	% 		Yhat, Z, Uhat, b, f, k, dfdz, dgdy, and dkdu.
	%
	% Outputs:
	%	v: The value of T at x given Tparams.

	% unpack Tparams
	ny = Tparams.ny;
	nz = Tparams.nz;
	nu = Tparams.nu;
	h = Tparams.h;
	yn = Tparams.yn;
	zn = Tparams.zn;
	Yhat = Tparams.Yhat;
	Z = Tparams.Z;
	Uhat = Tparams.Uhat;
	b = Tparams.b;
	f = Tparams.f;
	k = Tparams.k;
	dfdz = Tparams.dfdz;
	dgdy = Tparams.dgdy;
	dkdu = Tparams.dkdu;

	s = size(b,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= nz+nu
		error('Error (T1_lobattoIII): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (T1_lobattoIII): x has too many columns.')
	end

	% unpack x
	znp = x(1:nz);
	Unps = x(nz+1:nz+nu);

	Y0 = yn; % Y0 = Y1 = yn
	Z0 = Z(1:nz); % Z0 = Z1
	U0 = Uhat(1:nu);

	Y_arg = [Y0; Yhat(1:(s-1)*ny)]; % [Y0; Y2; ...; Ys]
	Z_arg = [Z0; Z(nz+1:nz*s)]; % [Z0; Z2; ...; Zs]
	U_arg = [U0; Uhat(2*nu+1:nu*s); Unps]; % [U0; U2; ...; U_{s-1}; U_s]

	ynp = Yhat((s-1)*ny+1:s*ny); 

	% Find T(x)
	v1 = znp - zn - h*kron(b',eye(nz)) * computeK(Y_arg,Z_arg,U_arg,k,ny,nz,nu,s); 
	v2 = (1/h) * dgdy(ynp) * f(ynp,znp); % Note the scaling by 1/h

	v = [v1; v2];
end


