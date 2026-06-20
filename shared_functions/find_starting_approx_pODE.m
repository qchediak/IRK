function [Ynp0 Znp0] = find_starting_approx_pODE(ynm, znm, yn, zn, Yn, Zn, starting_approx_str, SA_params)
	%function [Ynp0 Znp0] = find_starting_approx_pODE(ynm, znm, yn, zn, Yn, Zn, starting_approx_str, SA_params)
	%
	% Finds the starting approximations Ynp0 and Znp0 to Ynp and Znp
	% based % on the previous internal stages Yn and Zn and the past
	% step (ynm,znm) % (if starting_approx_str='good'), or chooses
	% trivial starting approximations 
	% (if starting_approx_str='trivial').
	%
	% Inputs:
	%	(ynm,znm): The past step y_{n-1}, z_{n-1}
	%	(yn,zn): The past step y_n, zn
	%	Yn, Zn: vector of the internal stages Y_{n,1},...,Y_{n,s} and
	%		Z_{n,1},...,Z_{n,s}.
	%	starting_approx_str: Indicates whether to use trivial or good
	%		starting approximations.
	%	SA_params: Structure array with fields: 
	%		h: stepsize of the first step (does NOT factor in r)
	%		(b,c,A,Ahat): PRK coefficients
	%		r: ratio of second stepsize / first stepsize
	% 		f,g: RHS of the partitioned ODE y'=f(y,z), z'=g(y,z)
	%		method_str: String describing the method.  
	%
	% Outputs:
	% 		Y0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)].
	% 		Z0 = [Z_{n+1,1}^(0); ...; Z_{n+1,s}^(0)].
	%
	% The above notation assumes that the method is going
	% from time t=t_n to t=t_{n+1}, but this need not be the case.

	% unpack params
	h = SA_params.h;
	b = SA_params.b;
	c = SA_params.c;
	A = SA_params.A;
	Ahat = SA_params.Ahat;
	r = SA_params.r;
	f = SA_params.f;
	g = SA_params.g;
	method_str = SA_params.method_str;

	ny = size(ynm,1);
	nz = size(znm,1);
	s = length(b);

	% Boolean value declaring if we have coefficients necessary to make
	% a good starting approximation.
	good_coefficients = false;

	% --------------------------
	% Find starred coefficients, if applicable
	% --------------------------

	if strcmpi(method_str, 'Lobatto IIIA') && s==2
		beta_star = (1/2) * [0 0; -2*r-r^2 r^2];
		%beta = (1/2) * [1 1; 1-r^2 1+2r+r^2];
		E = [0 1; 1 0];

		beta_hat_star = beta_star;
		good_coefficients = true;
	elseif strcmpi(method_str, 'Lobatto IIIB') && s==2
		beta_star = [-(1/2)*r 0; -(1/2)*r 0];
		E = [0 1; 1 0];

		beta_hat_star = beta_star;
		good_coefficients = true;
	elseif strcmpi(method_str, 'Gauss') && s==2
		beta_star11 = -(1/12)*r*(2*sqrt(3) + (-3+2*sqrt(3))*r);
		beta_star12 = (1/12)*(-3+2*sqrt(3))*r*(2+r);
		beta_star21 = -(1/12)*(3+2*sqrt(3))*r*(2+r);
		beta_star22 = (1/12)*r*(2*sqrt(3)+(3+2*sqrt(3))*r);
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];

		beta_hat_star = beta_star;
		good_coefficients = true;
	elseif strcmpi(method_str, 'Radau IIA') && s==2
		beta_star11 = (1/12) * (-r) * (r+4);
		beta_star12 = (1/12) * r^2;
		beta_star21 = (1/12) * (-3*r) * (3*r+4);
		beta_star22 = (1/12) * 9 * r^2;
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];

		beta_hat_star = beta_star;
		good_coefficients = true;
	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') 
		if s==2
			beta_star = [0 0; -r 0];
			beta_hat_star = [-(1/2)*r 0; -(1/2)*r 0];
			E = [0 1; 1 0];

			good_coefficients = true;
		elseif s==3
			%beta_star = [0 0 0;
			%	-r*c(2)-r^2*c(2)^2 r^2*c(2)^2 0;
			%	-r*c(3)-r^2*c(3)^2 r^2*c(3)^2 0];

			beta_star = (1/24) * [0 0 0;
				-r*(2*r^2+9*r+12) 4*r^2*(r+3) -r^2*(2*r+3);
				-4*r*(4*r^2+9*r+6) 16*r^2*(2*r+3) -4*r^2*(4*r+3)];

			beta_hat_star = (1/12) * [r^2*(r+3) -2*r^2*(r+2) r^2*(r+1);
				-2*r*(r^2+3*r+3) 4*r^2*(r+2) -2*r^2*(r+1);
				-r*(5*r^2+15*r+12) 10*r^2*(r+2) -5*r^2*(r+1)];

			E = [0 0 1; 0 1 0; 1 0 0];
			good_coefficients = true;
		end
	end

	% check that we have the necessary coefficients if we are to make a good starting approximation
	if ~good_coefficients && strcmpi(starting_approx_str, 'good')
		error('Error (find_starting_approx_pODE): Starting approximation coefficients are not available for this pair of method_str and s.')
	end

	if strcmpi(starting_approx_str, 'trivial')
		Ynp0 = kron(ones(s,1),yn); 
		Znp0 = kron(ones(s,1),zn); 
	elseif strcmpi(starting_approx_str, 'good') 
		% find unstarred coefficients
		beta = kron(ones(s,1), b') - beta_star*E;
		beta_hat = kron(ones(s,1), b') - beta_hat_star*E;

		% compute Ynp0 and Znp0
		% Note that computeF can be used for both partitioned ODEs and index 3 DAEs.
		Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,nz,s);
		Znp0 = kron(ones(s,1), znm) + h * kron(beta_hat, eye(nz)) * computeF(Yn,Zn,g,ny,nz,s);
	else
		error('Error (find_starting_approx_pODE): Invalid value of starting_approx_str or starting approximation coefficients not availble for this combination of method_str and s.')
	end
end


