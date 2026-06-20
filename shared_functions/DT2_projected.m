function J = DT2_projected(x,Tparams)
	% function J = DT2_projected(x,Tparams)
	% Jacobian of the function T2_projected, which scales by a factor of h.
	%
	% Note that the linear system changes with each step, because (yn,zn) are
	% 	updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unpsp]: The argument of T.
	% 	Tparams: Structured array with fields ny, nz, nu, h, zn, Ynp,
	% 		Znp, Unp, b, f, k, dfdz, dgdy, and dkdu.
	%
	% Outputs:
	%	J: The value of DT at x given Tparams.

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
		error('Error (DT1_lobattoIIIC): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DT1_lobattoIIIC): x has too many columns.')
	end

	% unpack x
	znp = x(1:nz);
	Unpsp = x(nz+1:nz+nu);

	Ynpsp = Ynp((s-1)*ny+1:s*ny); % Y_{n+1,s+1} = Y_{n+1,s}
	Znpsp = Znp((s-1)*nz+1:s*nz); % Z_{n+1,s+1} = Z_{n+1,s}
	ynp = Ynpsp; % y_{n+1} = Y_{n+1,s+1}.  For readability.

	% calculate the Jacobian as a 2x2 block matrix
	J11 = eye(nz);
	J12 = -h*b(s)*dkdu(Ynpsp,Znpsp,Unpsp);
	J21 = (1/h) * dgdy(ynp) * dfdz(ynp,znp); % scale by 1/h so the Jacobian remains invertible
	J22 = zeros(nu);

	J = [J11 J12; J21 J22];
end

