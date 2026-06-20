function v = T2_projected(x,Tparams)
	% function v = T2_projected(x,Tparams)
	% RHS of second linear system to solve for a projected method
	% applied to an index 3 DAE.  Does NOT scale the constraint equation
	% by a factor of 1/h.
	%
	% Note that the linear system changes with each step, because (yn,zn)
	% are updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unpsp]: The argument of T.
	% 	Tparams: Structured array with fields ny, nz, nu, h, zn, Ynp,
	% 		Znp, Unp, b, f, k, dfdz, dgdy, and dkdu.
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
	Ynp = Tparams.Ynp;
	Znp = Tparams.Znp;
	Unp = Tparams.Unp;
	b = Tparams.b;
	f = Tparams.f;
	k = Tparams.k;
	dfdz = Tparams.dfdz;
	dgdy = Tparams.dgdy;
	dkdu = Tparams.dkdu;

	s = size(b,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= nz+nu
		error('Error (T1_projected): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (T1_projected): x has too many columns.')
	end

	% unpack x
	znp = x(1:nz);
	Unpsp = x(nz+1:nz+nu);

	Ynpsp = Ynp((s-1)*ny+1:s*ny); % Y_{n+1,s+1} = Y_{n+1,s}
	Znpsp = Znp((s-1)*nz+1:s*nz); % Z_{n+1,s+1} = Z_{n+1,s}
	ynp = Ynpsp; % y_{n+1} = Y_{n+1,s+1}.  For readability.

	% Compute the vector [k(Y1,Z1,U1);...;k(Y_{s-1},Z_{s-1},U_{s-1})]
	% This avoids computing the ununsed value k(Y_s,Z_s,U_s)
	Y1 = Ynp(1:ny);
	Z1 = Znp(1:nz);
	U1 = Unp(1:nu);

	k_vec = k(Y1,Z1,U1);
	for i=2:s-1
		Y_temp = Ynp((i-1)*ny+1:i*ny);
		Z_temp = Znp((i-1)*nz+1:i*nz);
		U_temp = Unp((i-1)*nu+1:i*nu);

		k_vec = [k_vec; k(Y_temp,Z_temp,U_temp)];
	end

	% Find T(x)
	%v1 = znp - zn - h*kron([b(1:s-1)' 0],eye(nz)) * computeK(Ynp,Znp,Unp,k,ny,nz,nu,s) - h*b(s)*k(Ynpsp,Znpsp,Unpsp);
	v1 = znp - zn - h*kron(b(1:s-1)',eye(nz)) * k_vec - h*b(s)*k(Ynpsp,Znpsp,Unpsp);
	v2 = dgdy(ynp) * f(ynp,znp); 

	v = [v1; v2];
end


