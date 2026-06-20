function Ynp0 = find_starting_approx_ODE(ynm, yn, Yn, starting_approx_str, SA_params)
	%function Ynp0 = find_starting_approx_ODE(ynm, yn, Yn, starting_approx_str, SA_params)
	%
	% Finds the starting approximation Ynp0 to Ynp based
	% on the previous internal stages Yn and the past step ynm 
	% (if starting_approx_str='good'), or chooses trivial starting
	% approximations (if starting_approx_str='trivial').
	%
	% Inputs:
	%	ynm: The past step y_{n-1}
	%	yn: The past step y_n
	%	Yn: vector of the internal stages Y_{n,1},...,Y_{n,s}.
	%	starting_approx_str: Indicates whether to use trivial or good
	%		starting approximations.
	%	SA_params: Structure array with fields: 
	%		h: stepsize of the first step (does NOT factor in r)
	%		(b,c,A): IRK coefficients
	%		r: ratio of second stepsize / first stepsize
	% 		f: RHS of the ODE y'=f(y).
	%		method_str: String describing the method.  
	%
	% Outputs:
	% 		Y0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)].
	%
	% The above notation assumes that the method is going
	% from time t=t_n to t=t_{n+1}, but this need not be the case.

	% unpack params
	h = SA_params.h;
	b = SA_params.b;
	c = SA_params.c;
	A = SA_params.A;
	r = SA_params.r;
	f = SA_params.f;
	method_str = SA_params.method_str;

	ny = size(ynm,1);
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

		good_coefficients = true;
	elseif strcmpi(method_str, 'Lobatto IIIB') && s==2
		beta_star = [-(1/2)*r 0; -(1/2)*r 0];
		E = [0 1; 1 0];

		good_coefficients = true;
	elseif strcmpi(method_str, 'Gauss') && s==2
		beta_star11 = -(1/12)*r*(2*sqrt(3) + (-3+2*sqrt(3))*r);
		beta_star12 = (1/12)*(-3+2*sqrt(3))*r*(2+r);
		beta_star21 = -(1/12)*(3+2*sqrt(3))*r*(2+r);
		beta_star22 = (1/12)*r*(2*sqrt(3)+(3+2*sqrt(3))*r);
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];

		good_coefficients = true;
	elseif strcmpi(method_str, 'Radau IIA') && s==2
		beta_star11 = (1/12) * (-r) * (r+4);
		beta_star12 = (1/12) * r^2;
		beta_star21 = (1/12) * (-3*r) * (3*r+4);
		beta_star22 = (1/12) * 9 * r^2;
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];

		good_coefficients = true;
	end

	% check that we have the necessary coefficients if we are to make a good starting approximation
	if ~good_coefficients && strcmpi(starting_approx_str, 'good')
		error('Error (find_starting_approx_ODE): Starting approximation coefficients are not available for this pair of method_str and s.')
	end

	if strcmpi(starting_approx_str, 'trivial')
		Ynp0 = kron(ones(s,1),yn); 
	elseif strcmpi(starting_approx_str, 'good') 
		% find unstarred coefficients
		beta = kron(ones(s,1), b') - beta_star*E;

		% compute Ynp0
		Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF_ODE(Yn,f,ny,s);
	else
		error('Error (find_starting_approx_ODE): Invalid value of starting_approx_str or starting approximation coefficients not availble for this combination of method_str and s.')
	end
end


