function result = get_values_method(method_str, s)
	% function result = get_values_method(method_str, s)
	%
	% Returns some values (such as PRK coefficients) given a method_str.
	% Also performs a check on whether the combination of method_str
	% and s makes sense.  Works for ODEs, partitioned ODEs, and DAEs.
	%
	% Inputs:
	%	method_str: String describing the method.
	%	s: Number of stages.
	%
	% Outputs:
	%	result: Structured array with the fields:
	%		A_invertible: Boolean describing if A is invertible.
	%		Ahat_invertible: Boolean describing if Ahat is invertible.
	%		W: Inverse of A.  NaN if A is singular.
	%		What: Inverse of Ahat.  NaN if Ahat is singular.
	%		b: Column vector of weights.
	%		c: Column vector of nodes.
	%		A: PRK coefficient matrix.
	%		Ahat: PRK coefficient matrix.
	%		bstar: Weights of the reverse method.
	%		cstar: Nodes of the reverse method.
	%		Astar: PRK coefficient matrix of the reverse method.
	%		Ahatstar: PRK coefficient matrix of the reverse method.
	%		Wstar: Inverse of Astar.  NaN if Astar is singular.
	%		Whatstar: Inverse of Ahatstar.  NaN if Ahatstar is singular.
	%		reverse_method_str: String of the reverse method
	%			corresponding to method_str.
	%		number_nonlinear_systems: The number of nonlinear systems
	%			needed to solve the method over a single step.

	% ----------------------------------------------------------
	% check that the method and number of stages are compatible
	% ----------------------------------------------------------
	if strcmpi(method_str, 'Gauss')
		if ~(s==2 || s==3)
			error('Error (get_values_method): Gauss currently requires either s=2 or s=3')
		end
	elseif strcmpi(method_str, 'Radau IIA') || strcmpi(method_str, 'projected Radau IIA')
		if s==1 || s >= 4
			error('Error (get_values_method): Radau IIA currently requires either s=2 or s=3')
		end
	elseif strcmpi(method_str, 'Lobatto IIIC')
		if s==1 || s >= 4
			error('Error (get_values_method): Lobatto IIIC currently requires either s=2 or s=3')
		end
	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
		if s==1 || s >= 4
			error('Error (get_values_method): Lobatto IIIA-IIIB currently requires either s=2 or s=3')
		end
	elseif strcmpi(method_str, 'projected Lobatto IIIC')
		if s~=3
			error('Error (get_values_method): projected Lobatto IIIC currently requires s=3')
		end
	elseif strcmpi(method_str, 'projected Lobatto III')
		if s~=3
			error('Error (get_values_method): projected Lobatto III currently requires s=3')
		end
	end

	% ----------------------------------------------------------------
	% set the boolean values A_invertible and Ahat_invertible and the
	% inverse matrices 
	% ----------------------------------------------------------------
	if strcmpi(method_str, 'Gauss') 

		A_invertible = true;
		Ahat_invertible = true;

		Astar_invertible = true;
		Ahatstar_invertible = true;

		if s==2

			W11 = 3;
			W12 = -3+2*sqrt(3);
			W21 = -3-2*sqrt(3);
			W22 = 3;

			W = [W11 W12; W21 W22];
			What = W;

			Wstar = W;
			Whatstar = What;

		elseif s==3

			W11 = 5;
			W12 = -4+4*sqrt(5/3);
			W13 = 5-4*sqrt(5/3);

			W21 = -(5/6)*(3+sqrt(15));
			W22 = 2;
			W23 = (5/6)*(-3+sqrt(15));

			W31 = 5+4*sqrt(5/3);
			W32 = -(4/3)*(3+sqrt(15));
			W33 = 5;

			W = [W11 W12 W13; W21 W22 W23; W31 W32 W33];
			What = W;

			Wstar = W;
			Whatstar = What;
		else
			error('Error (get_values_method): Need to specify W=inverse(A) for Gauss')
		end

	elseif strcmpi(method_str, 'Radau IIA') || strcmpi(method_str, 'projected Radau IIA')

		A_invertible = true;
		Ahat_invertible = true;

		Astar_invertible = false;
		Ahatstar_invertible = false;

		if s==2
			W = [3/2 1/2; -9/2 5/2];
			What = W;
		elseif s==3
			W11 = 15*(4+sqrt(6));
			W12 = -36+29*sqrt(6);
			W13 = 12-8*sqrt(6);

			W21 = -36-29*sqrt(6);
			W22 = -15*(-4+sqrt(6));
			W23 = 12+8*sqrt(6);

			W31 = -30+80*sqrt(6);
			W32 = -10*(3+8*sqrt(6));
			W33 = 150;

			W = (1/30) * [W11 W12 W13; W21 W22 W23; W31 W32 W33];
			What = W;
		else
			error('Error (get_values_method): Need to specify W=inverse(A) for the Radau IIA method')
		end

	elseif strcmpi(method_str, 'Lobatto IIIC') || strcmpi(method_str, 'projected Lobatto IIIC')

		A_invertible = true;
		Ahat_invertible = true;

		Astar_invertible = false;
		Ahatstar_invertible = false;

		if s==2
			W = [1 1; -1 1];
			What = W;
		elseif s==3
			W = [3 4 -1; -1	0 1; 1 -4 3];
			What = W;
		else
			error('Error (get_values_method): Need to specify W=inverse(A) for the Lobatto IIIC method')
		end

	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') || strcmpi(method_str, 'Lobatto III') || strcmpi(method_str, 'projected Lobatto III') || strcmpi(method_str, 'Lobatto IIIA') || strcmpi(method_str, 'Lobatto IIIB')

		A_invertible = false;
		Ahat_invertible = false;

		Astar_invertible = false;
		Ahatstar_invertible = false;

	else

		error('Error (get_values_method): Invalid value of method_str.')

	end

	% ----------------------
	% set reverse_method_str
	% ----------------------
	if strcmpi(method_str, 'Gauss') || strcmpi(method_str, 'Lobatto IIIA-IIIB') || strcmpi(method_str, 'Lobatto IIIA') || strcmpi(method_str, 'Lobatto IIIB')
		reverse_method_str = method_str; % symmetric methods
	elseif strcmpi(method_str, 'Radau IIA')
		reverse_method_str = 'Radau I';
	elseif strcmpi(method_str, 'projected Radau IIA')
		reverse_method_str = 'projected Radau I';
	elseif strcmpi(method_str, 'Lobatto IIIC')
		reverse_method_str = 'Lobatto III';
	elseif strcmpi(method_str, 'projected Lobatto IIIC')
		reverse_method_str = 'projected Lobatto III';
	elseif strcmpi(method_str, 'Lobatto III')
		reverse_method_str = 'Lobatto IIIC';
	elseif strcmpi(method_str, 'projected Lobatto III')
		reverse_method_str = 'projected Lobatto IIIC';
	else
		error('Error (get_values_method): Invalid method_str.')
	end

	% ---------------------
	% set PRK coefficients
	% ---------------------
	if strcmpi(method_str, 'Gauss') 
		if s==2
			% 2-stage Gauss
			A = [1/4 1/4-sqrt(3)/6;
				1/4+sqrt(3)/6 1/4];
			b = [1/2; 1/2];
			c = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
			Ahat = A;
		elseif s==3
			% 3-stage Gauss
			A=[5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
				5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
				5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
			b = [5/18; 4/9; 5/18];
			c = [1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10];
			Ahat = A;
		end
	elseif strcmpi(method_str, 'Radau IIA') || strcmpi(method_str, 'projected Radau IIA')
		if s==2
			% 2-stage Radau IIA
			A = [5/12 -1/12;
				3/4 1/4];
			b = [3/4; 1/4];
			c = [1/3; 1];
		elseif s==3
			% 3-stage Radau IIA
			A11 =5*(88-7*sqrt(6));
			A12 = 296-169*sqrt(6);
			A13 = 8*(-2+3*sqrt(6));

			A21 = 296+169*sqrt(6);
			A22 = 5*(88+7*sqrt(6));
			A23 = -8*(2+3*sqrt(6));

			A31 = -50*(-16+sqrt(6));
			A32 = 50*(16+sqrt(6));
			A33 = 200;

			A = (1/1800) * [A11 A12 A13; A21 A22 A23; A31 A32 A33];
			b = [(16-sqrt(6))/36; (16+sqrt(6))/36; 1/9];
			c = [(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
		end
		Ahat = A;
	elseif strcmpi(method_str, 'Lobatto IIIA') 
		if s==2
			% 2-stage Lobatto IIIA
			A = [0 0; 1/2 1/2];
			Ahat = A;
			b = [1/2; 1/2];
			c = [0; 1];
		elseif s==3
			% 3-stage Lobatto IIIA
			A = [0 0 0;
				5/24 1/3 -1/24;
				1/6 2/3 1/6];
			Ahat = A;
			b = [1/6; 2/3; 1/6];
			c = [0; 1/2; 1];
		end
	elseif strcmpi(method_str, 'Lobatto IIIB')
		if s==2
			% 2-stage Lobatto IIIB
			A = [1/2 0; 1/2 0];
			Ahat = A;
			b = [1/2; 1/2];
			c = [0; 1];
		elseif s==3
			% 3-stage Lobatto IIIB
			A = [1/6 -1/6 0;
					1/6 1/3 0;
					1/6 5/6 0];
			Ahat = A;
			b = [1/6; 2/3; 1/6];
			c = [0; 1/2; 1];
		end
	elseif (strcmpi(method_str, 'Lobatto IIIC') || strcmpi(method_str, 'projected Lobatto IIIC')) 
		if s==2
			% 2-stage Lobatto IIIC
			A = [1/2 -1/2; 1/2 1/2];
			b = [1/2; 1/2];
			c = [0; 1];
			Ahat = A;
		elseif s==3
			% 3-stage Lobatto IIIC
			A = [1/6 -1/3 1/6;
				1/6 5/12 -1/12;
				1/6 2/3 1/6];
			b = [1/6; 2/3; 1/6];
			c = [0; 1/2; 1];
			Ahat = A;
		end
	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
		if s==2
			% 2-stage Lobatto IIIA-IIIB
			A = [0 0; 1/2 1/2];
			Ahat = [1/2 0; 1/2 0];
			b = [1/2; 1/2];
			c = [0; 1];
		elseif s==3
			% 3-stage Lobatto IIIA-IIIB
			A = [0 0 0;
				5/24 1/3 -1/24;
				1/6 2/3 1/6];

			Ahat = [1/6 -1/6 0;
					1/6 1/3 0;
					1/6 5/6 0];

			b = [1/6; 2/3; 1/6];
			c = [0; 1/2; 1];
		end
	elseif (strcmpi(method_str, 'Lobatto III') || strcmpi(method_str, 'projected Lobatto III')) && s==3
		% 3-stage Lobatto III
		A = [0 0 0; 
			1/4 1/4 0; 
			0 1 0];
		Ahat = A;
		b = [1/6; 2/3; 1/6];
		c = [0; 1/2; 1];
	else
		fprintf('Error (get_values_method): Invalid combination of method_str and s.')
	end

	% if a matrix is not invertible, set the "inverse" to nan
	if ~A_invertible
		W = nan;
	end

	if ~Ahat_invertible
		What = nan;
	end

	if ~Astar_invertible
		Wstar = nan;
	end

	if ~Ahatstar_invertible
		Whatstar = nan;
	end

	% Set the reverse coefficients
	E = flip(eye(s));
	bstar = E*b;
	cstar = E*(ones(s,1)-c);
	Astar = E*(kron(ones(s,1),b')-A)*E;
	Ahatstar = E*(kron(ones(s,1),b')-Ahat)*E;

	% For DAEs only:
	% Set the number of nonlinear systems based on the method.  
	% This makes the code easier to read.
	% Ignore this for ODEs and partitioned ODEs.
	if strcmpi(method_str, 'Lobatto III')
		number_nonlinear_systems = nan; % not necessary to set for DAEs.  Only a placeholder for ODEs.
	elseif strcmpi(method_str, 'Gauss') || strcmpi(method_str, 'Radau IIA') || strcmpi(method_str, 'Lobatto IIIC') || strcmpi(method_str, 'Lobatto IIIA') || strcmpi(method_str, 'Lobatto IIIB') 
		number_nonlinear_systems = 1;
	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') || strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Lobatto III') || strcmpi(method_str, 'projected Radau IIA')
		number_nonlinear_systems = 2;
	else
		error('Error (get_values_method): Invalid value of method_str')
	end

	result.A_invertible = A_invertible;
	result.Ahat_invertible = Ahat_invertible;
	result.b = b;
	result.c = c;
	result.A = A;
	result.Ahat = Ahat;
	result.W = W;
	result.What = What;

	result.Astar_invertible = Astar_invertible;
	result.Ahatstar_invertible = Ahatstar_invertible;
	result.bstar = bstar;
	result.cstar = cstar;
	result.Astar = Astar;
	result.Ahatstar = Ahatstar;
	result.Wstar = Wstar;
	result.Whatstar = Whatstar;

	result.number_nonlinear_systems = number_nonlinear_systems;
	result.reverse_method_str = reverse_method_str;

end



