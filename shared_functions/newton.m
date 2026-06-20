function [root iter] = newton(F,DF,x0,tol)
	% function [root iter] = newton(F,DF,x0,tol)
	% Solves F(x) = 0 by Newton's method.
	%
	% Inputs:
	%	F: RHS of the system 0=F(x) to solve.
	%	DF: The Jacobian of F.
	%	x0: The starting approximation.
	%	tol: tolerance.
	%
	% Outputs:
	%	root: The root x of F(x).
	%	iter: Number of Newton iterations.

	if size(x0,2) > 1
		error('Error (newton): x0 is not a column vector.')
	elseif size(F(x0),2) > 1
		error('Error (newton): F(x0) is not a column vector.')
	elseif size(DF(x0),1) ~= size(DF(x0),2)
		error('Error (newton): DF(x0) is not a square matrix.')
	elseif size(DF(x0),2) ~= size(F(x0),1)
		error('Error (newton): F(x0) and DF(x0) have incompatible dimensions.')
	end

	xn = x0;
	iter = 0;

	while max(abs(F(xn))) > tol
		iter = iter + 1;
		if iter > 40
			error('Newtons method did not converge')
		end

		delta = DF(xn) \ (-F(xn));
		xn = xn + delta;
	end

	if any(isnan(xn))
		error('Error (newton): The root contains one or more NaN values')
	else
		root = xn;
	end
end

