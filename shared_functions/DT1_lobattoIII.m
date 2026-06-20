function J = DT1_lobattoIII(x,Tparams)
	% function J = DT1_lobattoIII(x,Tparams)
	% Jacobian of the function T1.
	%
	% Note that the linear system changes with each step, because (yn,zn) are
	% 	updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unps]: The argument of T.
	% 	Tparams: Structured array with fields ny, nz, nu, h, yn, zn,
	% 		Yhat, Z, Uhat, b, f, k, dfdz, dgdy, and dkdu.
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
		error('Error (DT1_lobattoIII): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (DT1_lobattoIII): x has too many columns.')
	end

	% unpack x
	znp = x(1:nz);
	Unps = x(nz+1:nz+nu);

	Ynps = Yhat((s-2)*ny+1:(s-1)*ny); 
	Znps = Z((s-1)*nz+1:s*nz); 

	ynp = Yhat((s-1)*ny+1:s*ny); 

	% calculate the Jacobian as a 2x2 block matrix
	J11 = eye(nz);
	J12 = -h*b(s)*dkdu(Ynps,Znps,Unps);
	J21 = (1/h) * dgdy(ynp) * dfdz(ynp,znp); % scale by 1/h
	J22 = zeros(nu);

	J = [J11 J12; J21 J22];
end


