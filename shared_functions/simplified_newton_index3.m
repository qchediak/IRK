function [root iter] = simplified_newton_index3(F,x0,tol,params)
	% function [root iter] = simplified_newton_index3(F,x0,tol,params)
	%
	% A simplified Newton method with modified Newton iterations used
	% for solving the nonlinear system of an index 3 DAE.
	%
	% Importantly, instead of solving for \Delta U, this factors out an
	% h from the third column and solves for h \Delta U.
	%
	% Inputs:
	%	F: The nonlinear system to solve is F(x)=0.
	%	x0: The initial guess
	%	yn,zn,un: The previous step
	%	tol: Tolerance.
	%	params: Structured array consisting of derivatives of the
	%		functions, PRK coefficients, and dimensions.
	%
	% Outputs:
	%	root: The root of F(x)=0
	%	iter: The number of iterations to find the root.

	if size(x0,2) > 1
		error('Error (simplified_newton): x0 is not a column vector.')
	elseif size(F(x0),2) > 1
		error('Error (simplified_newton): F(x0) is not a column vector.')
	end

	% unpack params
	dfdz = params.dfdz;
	dkdu = params.dkdu;
	dgdy = params.dgdy;
	A = params.A;
	Ahat = params.Ahat;
	ny = params.ny;
	nz = params.nz;
	nu = params.nu;
	h = params.h;
	yn = params.yn;
	zn = params.zn;
	un = params.un;

	s = size(A,1);

	% make the approximate Jacobian matrix J 

	J11 = eye(ny*s);
	J12 = -h*kron(A,dfdz(yn,zn));
	J13 = zeros(ny*s,nu*s);

	J21 = zeros(s*nz,s*ny);
	J22 = eye(nz*s);
	J23 = -kron(Ahat,dkdu(yn,zn,un)); % the h has been factored out

	J31 = zeros(s*nu,s*ny);
	J32 = kron(A, dgdy(yn)*dfdz(yn,zn));
	J33 = zeros(nu*s);

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];

	xn = x0;
	iter = 0;

	while max(abs(F(xn))) > tol
		iter = iter + 1;
		if iter > 10
			error('simplified_newton did not converge')
		end

		delta = J \ (-F(xn));

		% Replace h \Delta U with \Delta U
		h_delta_U = delta(s*ny+s*nz+1:s*ny+s*nz+s*nu);
		delta_U = (1/h) * h_delta_U;
		delta(s*ny+s*nz+1:s*ny+s*nz+s*nu) = delta_U;

		xn = xn + delta;
	end

	if any(isnan(xn))
		error('Error (simplified_newton): The root contains one or more NaN values')
	else
		root = xn;
	end
end


