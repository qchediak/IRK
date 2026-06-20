function [root iter] = simplified_newton_index3_second_system(F,x0,tol,params)
	% function [root iter] = simplified_newton_index3_second_system(F,x0,tol,params)
	%
	% Solves F(x) = 0 by the simplified Newton method, using a
	% constant matrix J which approximates the Jacobian	and a modified
	% Newton iteration which factors out the stepsize h.
	%
	% Note that this assumes the constraint equation is NOT scaled by 1/h.
	%
	% Inputs:
	%	F: The nonlinear system to solve is F(x)=0.
	%	x0: The initial guess
	%	tol: Tolerance.
	%	params: Structured array with fields listed below.
	%
	% Outputs:
	%	root: The root of F(x)=0
	%	iter: The number of iterations to find the root.

	% unpack params
	ny = params.ny;
	nz = params.nz;
	nu = params.nu;
	h = params.h;
	yn = params.yn;
	zn = params.zn;
	Ynp = params.Ynp;
	Znp = params.Znp;
	Unp = params.Unp;
	b = params.b;
	f = params.f;
	k = params.k;
	dfdz = params.dfdz;
	dgdy = params.dgdy;
	dkdu = params.dkdu;

	s = size(b,1);

	if size(x0,2) > 1
		error('Error (simplified_newton): x0 is not a column vector.')
	elseif size(F(x0),2) > 1
		error('Error (simplified_newton): F(x0) is not a column vector.')
	end

	% unpack x0
	znp0 = x0(1:nz);
	Unpsp0 = x0(nz+1:nz+nu);

	Ynpsp = Ynp((s-1)*ny+1:s*ny); % Y_{n+1,s+1} = Y_{n+1,s}
	Znpsp = Znp((s-1)*nz+1:s*nz); % Z_{n+1,s+1} = Z_{n+1,s}
	ynp = Ynpsp; % y_{n+1} = Y_{n+1,s+1}.  For ease of reading.

	% make the approximate Jacobian J based on DT(x0)
	J11 = eye(nz);
	J12 = -b(s)*dkdu(Ynpsp,Znpsp,Unpsp0); % note the h has been factored out
	J21 = dgdy(ynp) * dfdz(ynp,znp0); 
	J22 = zeros(nu);

	J = [J11 J12; J21 J22];

	xn = x0;
	iter = 0;

	while max(abs(F(xn))) > tol
		iter = iter + 1;
		if iter > 10
			error('simplified_newton did not converge')
		end

		delta = J \ (-F(xn));

		% delta = [\Delta znp, h \Delta Usp].  Divide h \Delta Usp by h.
		h_delta_Usp = delta(nz+1:nz+nu);
		delta_Usp = (1/h) * h_delta_Usp;
		delta(nz+1:nz+nu) = delta_Usp;

		xn = xn + delta;
	end

	if any(isnan(xn))
		error('Error (simplified_newton_index3_second_system): The root contains one or more NaN values')
	else
		root = xn;
	end
end


