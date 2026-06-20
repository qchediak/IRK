function result = find_starting_approx(ynm, znm, yn, zn, un, Yn, Zn, Un, starting_approx_str, SA_params)
	% function result = find_starting_approx(ynm, znm, yn, zn, un, Yn, Zn, Un, starting_approx_str, SA_params)
	%
	% Finds the starting approximations Ynp0, Znp0, Unp0 to Ynp, Znp, Unp based
	% on the previous internal stages Yn, Zn, Un and the past step ynm, znm 
	% (if starting_approx_str='good'), or chooses trivial starting
	% approximations (if starting_approx_str='trivial').
	%
	% A third option, 'interpolation', gives starting approximations
	% based on the interpolating polynomial through the internal stages
	% (Y_i,Z_i,U_i), i=1,...,s.
	%
	% Inputs:
	%	ynm, znm: The past step y_{n-1}, z_{n-1}
	%	yn, zn, un: The past step y_n, z_n, u_n
	%	Yn: vector of the internal stages Y_{n,1},...,Y_{n,s}.
	%	Zn: vector of the internal stages Z_{n,1},...,Z_{n,s}.
	%	Un: vector of the internal stages U_{n,1},...,U_{n,s}.
	%	starting_approx_str: Indicates whether to use trivial or good
	%		starting approximations.
	%	SA_params: Structure array with fields: 
	%		h: stepsize of the first step (does NOT factor in r)
	%		(b,c,A,Ahat): PRK coefficients
	%		r: ratio of second stepsize / first stepsize
	% 		f: RHS of y'=f(y,z).
	% 		k: RHS of z'=f(y,z,u).
	%		nu: dimension of u
	%		method: String describing the method.  
	%		number_nonlinear_systems: The number of nonlinear systems.
	%
	% Outputs:
	%	result: A structure array.  If the number of nonlinear systems is 1, then 
	%		result has the fields Y0, Z0, U0, where
	% 			Y0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)]
	% 			Z0 = [Z_{n+1,1}^(0); ...; Z_{n+1,s}^(0)]
	% 			U0 = [U_{n+1,1}^(0); ...; U_{n+1,s}^(0)]
	%
	%		For the Lobatto IIIA-IIIB method, result has the fields
	%			Y0, Z0, U0, znp0, Unps0, where
	%			Y0 = [Y_{n+1,2}^{(0)}; ...; Y_{n+1,s}^{(0)}]
	%			Z0 = [Z_{n+1,1}^{(0)}; ...; Z_{n+1,s}^{(0)}]
	%			U0 = [U_{n+1,1}^{(0)}; ...; U_{n+1,s-1}^{(0)}]
	%			z0 = z_{n+1}^{(0)}
	%			Us0 = U_{n+1,s}^{(0)}
	%
	%		For the projected Lobatto IIIC or projected Radau IIA methods, 
	%			result has the fields Y0, Z0, U0, znp0, Unpsp0, where
	% 			Y0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)]
	% 			Z0 = [Z_{n+1,1}^(0); ...; Z_{n+1,s}^(0)]
	% 			U0 = [U_{n+1,1}^(0); ...; U_{n+1,s}^(0)]
	%			z0 = z_{n+1}^{(0)}
	%			Usp0 = U_{n+1,s+1}^{(0)}
	%
	%		For the projected Lobatto III method, 
	%			result has the fields Yhat0, Z0, Uhat0, znp0, Us0, where
	% 			Yhat0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)]
	% 			Z0 = [Z_{n+1,1}^(0); ...; Z_{n+1,s}^(0)]
	% 			Uhat0 = [U_{n+1,1}^(0); ...; U_{n+1,s}^(0)]
	%			znp0 = z_{n+1}^{(0)}
	%			Us0 = U_{n+1,s+1}^{(0)}
	%
	%	The above notation assumes that the method is going
	%	from time t=t_n to t=t_{n+1}, but this need not be the case.

	% unpack params
	h = SA_params.h;
	b = SA_params.b;
	c = SA_params.c;
	A = SA_params.A;
	Ahat = SA_params.Ahat;
	r = SA_params.r;
	f = SA_params.f;
	k = SA_params.k;
	nu = SA_params.nu;
	method_str = SA_params.method_str;
	number_nonlinear_systems = SA_params.number_nonlinear_systems;

	ny = size(ynm,1);
	nz = size(znm,1);
	s = length(b);

	% --------------------------
	% Find starred coefficients
	% --------------------------

	% Check that ynm, znm, unm are not NaN with good starting approximations
	if anynan([ynm; znm]) && strcmpi(starting_approx_str, 'good')
		error('Error (find_starting_approx): Cannot use good starting approximations if any of the values in ynm, znm, or unm is NaN.')
	end

	% projected Lobatto III has no preset good starting approximation
	if strcmpi(method_str, 'projected Lobatto III') && strcmpi(starting_approx_str, 'good')
		error('Error (find_starting_approx): This method has no preset good starting approximation')
	end

	if strcmpi(starting_approx_str, 'interpolation')

		% This currently requires that s=3 and tnm=0; i.e. the
		% interpolating polynomial interpolates over the time interval
		% [tnm,tn]=[0,h] and is extrapolated to the time interval
		% [tn,tnp]=[h,(1+r)h].
		if s ~= 3
			error('Error (find_starting_approx): The starting approximations based on interpolation currently only work for s=3')
		end

		tnm = 0; 

		Yn1 = Yn(1:ny);
		Yn2 = Yn(ny+1:2*ny);
		Yn3 = Yn(2*ny+1:3*ny);

		Zn1 = Zn(1:nz);
		Zn2 = Zn(nz+1:2*nz);
		Zn3 = Zn(2*nz+1:3*nz);

		Un1 = Un(1:nu);
		Un2 = Un(nu+1:2*nu);
		Un3 = Un(2*nu+1:3*nu);

		tau = @(t) (1/h) * (t-tnm);

		% Lagrange basis polynomials
		L1 = @(T) ( 1 / ( (c(1)-c(2))*(c(1)-c(3)) ) ) * (T-c(2))*(T-c(3));
		L2 = @(T) ( 1 / ( (c(2)-c(1))*(c(2)-c(3)) ) ) * (T-c(1))*(T-c(3));
		L3 = @(T) ( 1 / ( (c(3)-c(1))*(c(3)-c(2)) ) ) * (T-c(1))*(T-c(2));

		% interpolating polynomials 
		rhoY = @(t) L1(tau(t))*Yn1 + L2(tau(t))*Yn2 + L3(tau(t))*Yn3;
		rhoZ = @(t) L1(tau(t))*Zn1 + L2(tau(t))*Zn2 + L3(tau(t))*Zn3;
		rhoU = @(t) L1(tau(t))*Un1 + L2(tau(t))*Un2 + L3(tau(t))*Un3;

	elseif strcmpi(starting_approx_str, 'good')

		if strcmpi(method_str, 'Gauss') && s==3

			beta_star11 = (1/180)*r*(-30 + 3*(-25 + 6*sqrt(15))*r + 2*(-35 + 9*sqrt(15))*r^2);
			beta_star12 = -(1/45)*r*(3*(-5 + sqrt(15)) + 15*(-4 + sqrt(15))*r + (-35 + 9*sqrt(15))*r^2);
			beta_star13 = (1/180)*r*(30*(-4 + sqrt(15)) + 3*(-55 + 14*sqrt(15))*r + 2*(-35 + 9*sqrt(15))*r^2);

			beta_star21 = -(1/72)*r*(6*(5 + sqrt(15)) + 3*(10 + sqrt(15))*r + 10*r^2);
			beta_star22 = (1/18)*r*(6 + 15*r + 5*r^2);
			beta_star23 = (1/72)*r*(6*(-5 + sqrt(15)) + 3*(-10 + sqrt(15))*r - 10*r^2);

			beta_star31 = -(1/180)*r*(30*(4 + sqrt(15)) + 3*(55 + 14*sqrt(15))*r + 2*(35 + 9*sqrt(15))*r^2);
			beta_star32 = (1/45)*r*(3*(5 + sqrt(15)) + 15*(4 + sqrt(15))*r + (35 + 9*sqrt(15))*r^2);
			beta_star33 = -(1/180)*r*(30 + 3*(25 + 6*sqrt(15))*r + 2*(35 + 9*sqrt(15))*r^2);

			beta_star = [beta_star11 beta_star12 beta_star13; 
				beta_star21 beta_star22 beta_star23;
				beta_star31 beta_star32 beta_star33];

			E = [0 0 1; 0 1 0; 1 0 0];

			beta = kron(ones(s,1), b') - beta_star*E;
	
			beta_hat = beta;

			D_star11 = 5+sqrt(15)-(-7+sqrt(15))*r-2*(-4+sqrt(15))*r^2;
			D_star12 = 4*(1+r)*(-1+(-4+sqrt(15))*r);
			D_star13 = -((1+r)*(-5+sqrt(15)+2*(-4+sqrt(15))*r));

			D_star21 = (1+r)*(5+sqrt(15)+5*r);
			D_star22 = -2*(2+10*r+5*r^2);
			D_star23 = -((-5+sqrt(15)-5*r)*(1+r));

			D_star31 = (1+r)*(5+sqrt(15)+2*(4+sqrt(15))*r);
			D_star32 = -4*(1+r)*(1+(4+sqrt(15))*r);
			D_star33 = 5-sqrt(15)+(7+sqrt(15))*r+2*(4+sqrt(15))*r^2;

			D_star = (1/6) * [D_star11 D_star12 D_star13; 
				D_star21 D_star22 D_star23; 
				D_star31 D_star32 D_star33];

			D = D_star*E;

%			D_star11 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(1);
%			D_star12 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(1);
%			D_star13 = 0;
%			D_star21 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(2);
%			D_star22 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(2);
%			D_star23 = 0;
%			D_star31 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(3);
%			D_star32 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(3);
%			D_star33 = 0;
%			D_star = [D_star11 D_star12 D_star13; 
%				D_star21 D_star22 D_star23; 
%				D_star31 D_star32 D_star33];

		elseif strcmpi(method_str, 'Gauss') && s==2

			beta_star11 = -(1/12)*r*(2*sqrt(3) + (-3+2*sqrt(3))*r);
			beta_star12 = (1/12)*(-3+2*sqrt(3))*r*(2+r);
			beta_star21 = -(1/12)*(3+2*sqrt(3))*r*(2+r);
			beta_star22 = (1/12)*r*(2*sqrt(3)+(3+2*sqrt(3))*r);
			beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
			E = [0 1; 1 0];
			beta = kron(ones(s,1), b') - beta_star*E;

			beta_hat11 = 1+r*c(1);
			beta_hat12 = 0;
			beta_hat21 = 1+r*c(2);
			beta_hat22 = 0;
			beta_hat = [beta_hat11 beta_hat12; beta_hat21 beta_hat22];

			D11 = 1;
			D12 = 0;
			D21 = 1;
			D22 = 0;
			D = [D11 D12; D21 D22];

		elseif strcmpi(method_str, 'Radau IIA') && s==2

			beta_star11 = (1/12) * (-r) * (r+4);
			beta_star12 = (1/12) * r^2;
			beta_star21 = (1/12) * (-3*r) * (3*r+4);
			beta_star22 = (1/12) * 9 * r^2;
			beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
			E = [0 1; 1 0];
			beta = kron(ones(s,1), b') - beta_star*E;

			beta_hat11 = 1+r*c(1);
			beta_hat12 = 0;
			beta_hat21 = 1+r*c(2);
			beta_hat22 = 0;
			beta_hat = [beta_hat11 beta_hat12; beta_hat21 beta_hat22];

			D11 = 1;
			D12 = 0;
			D21 = 1;
			D22 = 0;
			D = [D11 D12; D21 D22];

		elseif (strcmpi(method_str, 'Lobatto IIIC') || strcmpi(method_str, 'projected Lobatto IIIC')) && s==3

			% shared first nonlinear system
			beta_star11 = -4*r^3;
			beta_star12 = 8*r^3;
			beta_star13 = -4*r^3;
			beta_star21 = -r*(r^2+9*r+12);
			beta_star22 = 2*r^2*(r+6);
			beta_star23 = -r^2*(r+3);
			beta_star31 = -4*r*(4*r^2+9*r+6);
			beta_star32 = 16*r^2*(2*r+3);
			beta_star33 = -4*r^2*(4*r+3);
			beta_star = (1/24) * [beta_star11 beta_star12 beta_star13; 
				beta_star21 beta_star22 beta_star23; 
				beta_star31 beta_star32 beta_star33];

			E = [0 0 1; 0 1 0; 1 0 0];

			beta = kron(ones(s,1), b') - beta_star*E;

			theta_star11 = -2*r^2;
			theta_star12 = 4*(r^2-1);
			theta_star13 = 2-2*r^2;
			theta_star21 = r^2;
			theta_star22 = -5*r^2 - 12*r - 4;
			theta_star23 = 4*r^2 + 6*r + 2;
			theta_star31 = -2*r^2;
			theta_star32 = -4*(2*r^2 + 6*r + 1);
			theta_star33 = 2*(5*r^2 + 6*r + 1);

			theta_star = (1/12) * [theta_star11 theta_star12 theta_star13;
									theta_star21 theta_star22 theta_star23;
									theta_star31 theta_star32 theta_star33];

			e1 = eye(1,s);

			beta_hat_star = theta_star + b(s)*kron(ones(s,1),e1);
			beta_hat = kron(ones(s,1), b') - beta_hat_star*E;

			D_star11 = 0;
			D_star12 = 2;
			D_star13 = -1;

			D_star21 = 0;
			D_star22 = r+2;
			D_star23 = -r-1;

			D_star31 = r;
			D_star32 = 2;
			D_star33 = -r-1;

			D_star = [D_star11 D_star12 D_star13; 
				D_star21 D_star22 D_star23;
				D_star31 D_star32 D_star33];

			D = D_star*E;

			% second nonlinear system for the projected method only
			if strcmpi(method_str, 'projected Lobatto IIIC')
				theta_star41 = 0;
				theta_star42 = -r^2 - 2*r - 1/3;
				theta_star43 = r^2 + r + 1/6;

				theta_star4 = [theta_star41 theta_star42 theta_star43];
				beta_hat_star4 = theta_star4 + b(s)*e1;

				beta_hat4 = b' - beta_hat_star4*E;

				D_star41 = 0;
				D_star42 = 2*r+2;
				D_star43 = -2*r-1;

				D_star4 = [D_star41 D_star42 D_star43];

				D4 = D_star4*E;
			end

		elseif strcmpi(method_str, 'projected Radau IIA') && s==3

			e1 = eye(1,s);
			E = [0 0 1; 0 1 0; 1 0 0];

			beta_star11 = 20*r*(-45*(-3+sqrt(6))+r*(162-63*sqrt(6)+(57-23*sqrt(6))*r));
			beta_star12 = 5*r^2*(117-63*sqrt(6)+(162-68*sqrt(6))*r);
			beta_star13 = 25*r^2*(-153+63*sqrt(6)-78*r+32*sqrt(6)*r);

			beta_star21 = 4*r*(45*(9+sqrt(6))+r*(378+117*sqrt(6)+123*r+47*sqrt(6)*r));
			beta_star22 = -5*r^2*(45*(7+3*sqrt(6))+2*(81+34*sqrt(6))*r);
			beta_star23 = r^2*(63+207*sqrt(6)+318*r+152*sqrt(6)*r);

			beta_star31 = -100*(-6+sqrt(6))*r*(9+2*r*(9+5*r));
			beta_star32 = -250*r^2*(9+9*sqrt(6)+10*sqrt(6)*r);
			beta_star33 = 50*r^2*(-171+81*sqrt(6)-120*r+70*sqrt(6)*r);

			beta_star = (1 / (900 * (-6 + sqrt(6)))) * [beta_star11 beta_star12 beta_star13; 
					beta_star21 beta_star22 beta_star23; 
					beta_star31 beta_star32 beta_star33];

			beta = kron(ones(s,1), b') - beta_star*E;

			%theta_star10 = (100*(-6723+1822*sqrt(6));
			theta_star11 = 294800-53450*sqrt(6)+r*(188160-24490*sqrt(6)+(-554+6731*sqrt(6))*r);
			theta_star12 = 625*(-316+199*sqrt(6))+250*(-8406+3959*sqrt(6))*r+(-2131946+879519*sqrt(6))*r^2;
			theta_star13 = -625*(-920+405*sqrt(6)+2*r*(-4254+1781*sqrt(6)+(-1706+709*sqrt(6))*r));

			%theta_star20 = 100*(-6723+1822*sqrt(6));
			theta_star21 = 25*(11792-2138*sqrt(6)+r*(11856+3866*sqrt(6)+(5494+2259*sqrt(6))*r));
			theta_star22 = 25*(-7900+4975*sqrt(6)+r*(5100+19850*sqrt(6)+2806*r+5191*sqrt(6)*r));
			theta_star23 = -625*(-920+405*sqrt(6)+2*r*(-810+515*sqrt(6)+(166+149*sqrt(6))*r));

			%theta_star30 = 100*(-6723+1822*sqrt(6));
			theta_star31 = 50*(1+r)*(5896-1069*sqrt(6)+(6218+2873*sqrt(6))*r);
			theta_star32 = 25*(-7900+4975*sqrt(6)+2*r*(50*(-987+743*sqrt(6))+(-93718+45877*sqrt(6))*r));
			theta_star33 = -3125*(-184+81*sqrt(6)+4*r*(-633+287*sqrt(6)+5*(-70+39*sqrt(6))*r));

			theta_star = ( 1 / (900*(-6723+1822*sqrt(6))) ) * [theta_star11 theta_star12 theta_star13;
									theta_star21 theta_star22 theta_star23;
									theta_star31 theta_star32 theta_star33];
			
			beta_hat_star = theta_star + b(s)*kron(ones(s,1),e1);
			beta_hat = kron(ones(s,1), b') - beta_hat_star*E;

			D_star11 = -5*(107142+27743*sqrt(6)+(19076+9985*sqrt(6))*r);
			D_star12 = 5*(25*(5472-9437*sqrt(6))+(646228-288837*sqrt(6))*r);
			D_star13 = 5*(-757740+385015*sqrt(6)+(-627152+298822*sqrt(6))*r);

			D_star21 = -535710-138715*sqrt(6)-449476*r-186139*sqrt(6)*r;
			D_star22 = 125*(5472-9437*sqrt(6))+(176420-592295*sqrt(6))*r;
			D_star23 = -3788700+1925075*sqrt(6)+273056*r+778434*sqrt(6)*r;

			D_star31 = -5*(107142+27743*sqrt(6)+2*(68107+29508*sqrt(6))*r);	
			D_star32 = -25*(5*(-5472+9437*sqrt(6))+2*(-85189+50912*sqrt(6))*r);
			D_star33 = 5*(-757740+385015*sqrt(6)+(-715676+568136*sqrt(6))*r);

			D_star = (1 / (5*121347*(-6+sqrt(6)))) * [D_star11 D_star12 D_star13; 
				D_star21 D_star22 D_star23;
				D_star31 D_star32 D_star33];

			D = D_star*E;

			% second nonlinear system
			theta_star41 = ((1+r)*(5896-1069*sqrt(6)+(6218+2873*sqrt(6))*r))/(18*(-6723+1822*sqrt(6)));
			theta_star42 = (-7900+4975*sqrt(6)+100*(-987+743*sqrt(6))*r+2*(-93718+45877*sqrt(6))*r^2)/(36*(-6723+1822*sqrt(6)));
			theta_star43 = -((125*(-309+151*sqrt(6)+6*(-692+363*sqrt(6))*r+20*(-93+82*sqrt(6))*r^2))/(54*(-9802+1403*sqrt(6)))); 

			theta_star4 = [theta_star41 theta_star42 theta_star43];

			beta_hat_star4 = theta_star4 + b(s)*e1;
			beta_hat4 = b' - beta_hat_star4*E;

			D_star41 = (6057+902*sqrt(6)+6218*r+2873*sqrt(6)*r)/(60507-16398*sqrt(6));
			D_star42 = (25*(987-743*sqrt(6))+(93718-45877*sqrt(6))*r)/(9*(-6723+1822*sqrt(6))); 
			D_star43 = (125*(3*(-692+363*sqrt(6))+20*(-93+82*sqrt(6))*r))/(27*(-9802+1403*sqrt(6)));

			D_star4 = [D_star41 D_star42 D_star43];

			D4 = D_star4*E;

		elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') && s==2

			beta11_star = 0; % only needed for matrix calculation
			beta12_star = 0; % only needed for matrix calculation
			beta21_star = -r;
			beta22_star = 0;

			beta11_hat_star = -(1/2) * r;
			beta12_hat_star = 0;
			beta21_hat_star = -(1/2) * r;
			beta22_hat_star = 0;
			beta31_hat_star = -r;
			beta32_hat_star = 0;

			D11_star = 1/2;
			D12_star = 1/2;
			D21_star = 1/2;
			D22_star = 1/2;

			beta_star = [beta11_star beta12_star; 
							beta21_star beta22_star];

			beta_hat_star = [beta11_hat_star beta12_hat_star; 
								beta21_hat_star beta22_hat_star;
								beta31_hat_star beta32_hat_star];

			D_star = [D11_star D12_star; 
						D21_star D22_star];

			E = [0 1; 1 0];
			beta = kron(ones(s,1), b') - beta_star*E;
			beta_hat = kron(ones(s+1,1), b') - beta_hat_star*E;
			D = D_star*E;

		elseif strcmpi(method_str, 'Lobatto IIIA-IIIB') && s==3

			beta_star = [0 0 0;
				-r*c(2)*(1+r*c(2)) r^2*c(2)^2 0;
				-r*c(3)*(1+r*c(3)) r^2*c(3)^2 0];

			beta_hat_star = [-r*c(1)-r^2*Ahat(1,2) r^2*Ahat(1,2) 0;
						-r*c(2)-r^2*Ahat(2,2) r^2*Ahat(2,2) 0;
						-r*c(3)-r^2*Ahat(3,2) r^2*Ahat(3,2) 0;
						-r^2-r r^2 0];

			D_star = [1+2*r*c(1) -2*r*c(1) 0;
				1+2*r*c(2) -2*r*c(2) 0;
				1+2*r*c(3) -2*r*c(3) 0];

			E = [0 0 1; 
				0 1 0; 
				1 0 0];

			beta = kron(ones(s,1),b') - beta_star*E;
			beta_hat = kron(ones(s+1,1),b') - beta_hat_star*E;
			D = D_star*E;

		else

			error('Error (find_SA): Invalid method_str and s value combination')

		end
	end

	% --------------------------------
	% Compute starting approximations
	% --------------------------------

	if number_nonlinear_systems==1

		if strcmpi(starting_approx_str, 'trivial')
			Ynp0 = kron(ones(s,1),yn); 
			Znp0 = kron(ones(s,1),zn); 
			Unp0 = kron(ones(s,1),un);
		elseif strcmpi(starting_approx_str, 'good')
			Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,nz,s);
			Znp0 = kron(ones(s,1), znm) + h * kron(beta_hat, eye(nz)) * computeK(Yn,Zn,Un,k,ny,nz,nu,s);
			Unp0 = kron(D, eye(nu)) * Un;
		elseif strcmpi(starting_approx_str, 'interpolation')

			% Recall this currently assumes s=3 and tnm=0
			Ynp10 = rhoY(tnm + (1+r*c(1))*h);
			Ynp20 = rhoY(tnm + (1+r*c(2))*h);
			Ynp30 = rhoY(tnm + (1+r*c(3))*h);

			Znp10 = rhoZ(tnm + (1+r*c(1))*h);
			Znp20 = rhoZ(tnm + (1+r*c(2))*h);
			Znp30 = rhoZ(tnm + (1+r*c(3))*h);

			Unp10 = rhoU(tnm + (1+r*c(1))*h);
			Unp20 = rhoU(tnm + (1+r*c(2))*h);
			Unp30 = rhoU(tnm + (1+r*c(3))*h);

			Ynp0 = [Ynp10; Ynp20; Ynp30];
			Znp0 = [Znp10; Znp20; Znp30];
			Unp0 = [Unp10; Unp20; Unp30];

		else
			error('Error (find_starting_approx): Invalid value of starting_approx_str')
		end

		result.Y0 = Ynp0;
		result.Z0 = Znp0;
		result.U0 = Unp0;

	elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')

		if strcmpi(starting_approx_str, 'trivial')
			Ynp0 = kron(ones(s-1,1),yn); 
			Znp0 = kron(ones(s,1),zn); 
			Unp0 = kron(ones(s-1,1),un);
			znp0 = zn;
			Unps0 = un;
		elseif strcmpi(starting_approx_str, 'good')
			Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,nz,s);
			Znp0 = kron(ones(s+1,1), znm) + h * kron(beta_hat, eye(nz)) * computeK(Yn,Zn,Un,k,ny,nz,nu,s);
			Unp0 = kron(D, eye(nu)) * Un;

			% Find znp0 and Unps0 from Znp0 and Unp0
			znp0 = Znp0(s*nz+1:(s+1)*nz);
			Unps0 = Unp0((s-1)*nu+1:s*nu);

			% Keep only the relevant parts of Ynp0, Znp0, and Unp0
			Ynp0 = Ynp0(ny+1:ny*s); 	% [Y_2^0; ...; Y_s^0]
			Znp0 = Znp0(1:s*nz); 		% [Z_1^0; ...; Z_s^0]
			Unp0 = Unp0(1:(s-1)*nu); 	% [U_1^0; ...; U_{s-1}^0]

		elseif strcmpi(starting_approx_str, 'interpolation')

			% Recall this currently assumes s=3 and tnm=0
			Ynp20 = rhoY(tnm + (1+r*c(2))*h);
			Ynp30 = rhoY(tnm + (1+r*c(3))*h);

			Znp10 = rhoZ(tnm + (1+r*c(1))*h);
			Znp20 = rhoZ(tnm + (1+r*c(2))*h);
			Znp30 = rhoZ(tnm + (1+r*c(3))*h);

			Unp10 = rhoU(tnm + (1+r*c(1))*h);
			Unp20 = rhoU(tnm + (1+r*c(2))*h);
			Unp30 = rhoU(tnm + (1+r*c(3))*h);

			Ynp0 = [Ynp20; Ynp30];
			Znp0 = [Znp10; Znp20; Znp30];
			Unp0 = [Unp10; Unp20];
			znp0 = rhoZ(tnm + (1+r)*h);
			Unps0 = Unp30;

		else
			error('Error (find_starting_approx): Invalid value of starting_approx_str')
		end

		result.Y0 = Ynp0;
		result.Z0 = Znp0;
		result.U0 = Unp0;
		result.z0 = znp0;
		result.Us0 = Unps0;

	elseif strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')

		if strcmpi(starting_approx_str, 'trivial')
			Ynp0 = kron(ones(s,1),yn); 
			Znp0 = kron(ones(s,1),zn); 
			Unp0 = kron(ones(s,1),un);
			znp0 = zn;
			Unpsp0 = un;
		elseif strcmpi(starting_approx_str, 'good')
			% Compute Ynp0, Znp0, Unp0
			Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,nz,s);
			Znp0 = kron(ones(s,1), znm) + h * kron(beta_hat, eye(nz)) * computeK(Yn,Zn,Un,k,ny,nz,nu,s);
			Unp0 = kron(D, eye(nu)) * Un;

			% Compute znp0 and Unpsp0
			znp0 = znm + h * kron(beta_hat4, eye(nz)) * computeK(Yn,Zn,Un,k,ny,nz,nu,s);
			Unpsp0 = kron(D4, eye(nu)) * Un;

		elseif strcmpi(starting_approx_str, 'interpolation')

			% Recall this currently assumes s=3 and tnm=0
			Ynp10 = rhoY(tnm + (1+r*c(1))*h);
			Ynp20 = rhoY(tnm + (1+r*c(2))*h);
			Ynp30 = rhoY(tnm + (1+r*c(3))*h);

			Znp10 = rhoZ(tnm + (1+r*c(1))*h);
			Znp20 = rhoZ(tnm + (1+r*c(2))*h);
			Znp30 = rhoZ(tnm + (1+r*c(3))*h);

			Unp10 = rhoU(tnm + (1+r*c(1))*h);
			Unp20 = rhoU(tnm + (1+r*c(2))*h);
			Unp30 = rhoU(tnm + (1+r*c(3))*h);

			Ynp0 = [Ynp10; Ynp20; Ynp30];
			Znp0 = [Znp10; Znp20; Znp30];
			Unp0 = [Unp10; Unp20; Unp30];
			znp0 = rhoZ(tnm + (1+r)*h);
			Unpsp0 = Unp30;

		else
			error('Error (find_starting_approx): Invalid value of starting_approx_str')
		end

		result.Y0 = Ynp0;
		result.Z0 = Znp0;
		result.U0 = Unp0;
		result.z0 = znp0;
		result.Usp0 = Unpsp0;

	elseif strcmpi(method_str, 'projected Lobatto III')

		if strcmpi(starting_approx_str, 'trivial')
			Yhat0 = kron(ones(s,1),yn); 
			Z0 = kron(ones(s,1),zn); 
			Uhat0 = kron(ones(s,1),un);
			znp0 = zn;
			Us0 = un;

		elseif strcmpi(starting_approx_str, 'interpolation')

			% Recall this currently assumes s=3 and tnm=0
			Ynp20 = rhoY(tnm + (1+r*c(2))*h);
			Ynp30 = rhoY(tnm + (1+r*c(3))*h);

			Znp10 = rhoZ(tnm + (1+r*c(1))*h);
			Znp20 = rhoZ(tnm + (1+r*c(2))*h);
			Znp30 = rhoZ(tnm + (1+r*c(3))*h);

			Unp10 = rhoU(tnm + (1+r*c(1))*h);
			Unp20 = rhoU(tnm + (1+r*c(2))*h);
			Unp30 = rhoU(tnm + (1+r*c(3))*h);

			ynp0 = Ynp30;
			Unp00 = Unp10;

			Yhat0 = [Ynp20; Ynp30; ynp0];
			Z0 = [Znp10; Znp20; Znp30];
			Uhat0 = [Unp00; Unp10; Unp20];
			znp0 = Znp30;
			Us0 = Unp30;

		else
			error('Error (find_starting_approx): Projected Lobatto III has no preset good starting approximation')
		end

		result.Yhat0 = Yhat0;
		result.Z0 = Z0;
		result.Uhat0 = Uhat0;
		result.znp0 = znp0;
		result.Us0 = Us0;

	else

		error('Error (find_starting_approx: Invalid value of number_nonlinear_systems or method_str')

	end
end



