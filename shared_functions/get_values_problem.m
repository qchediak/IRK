function result = get_values_problem(problem_str)
	% function result = get_values_problem(problem_str)
	%
	% Returns values based on problem_str.  For use when solving the index 3 DAE
	% y'=f(y,z), z'=k(y,z,u), 0=g(y).
	%
	% Inputs: 
	%	problem_str: String describing the index 3 DAE to be solved.
	%
	% Outputs:
	%	result: Structured array with the following fields:
	%		y0, z0, u0: Initial values
	%		f, k, g, dfdy, dfdz, dkdy, dkdz, dkdu, dgdy: Functions from
	%			the index 3 DAE and their derivatives.
	%		true_solution: A structured array with the following fields:
	%			y_true, z_true, u_true: True solution for y, z, and u.
	%				Returns NaN if unknown.

	% To try a custom problem, define (y0,z0,u0) and f, k, g, and their
	% derivatives in a separate branch below.  If the true solution
	% is known, it could also be defined here.

	% Note that if the iterative scheme is simplified Newton 3, then the
	% only derivatives that are needed are dfdz, dkdu, and dgdy.

	if strcmpi(problem_str, 'pendulum')

		% parameters related to the pendulum
		g0 = 9.8; % Graviational constant.
		alpha0 = [0;g0];
		l = 1; % length of rod

		% initial values.  
		y0 = [3/5; -4/5]; 
		z0 = [-4/2; -3/2];
		u0 = (1/l^2) * (dot(z0,z0) - dot(y0,alpha0)); % satisfies hidden constraint.  

		f = @(y,z) z;
		k = @(y,z,u) -u*y - alpha0;
		g = @(y) dot(y,y) - l^2;

		dfdy = @(y,z) zeros(2);
		dfdz = @(y,z) eye(2);
		dkdy = @(y,z,u) -u*eye(2);
		dkdz = @(y,z,u) zeros(2);
		dkdu = @(y,z,u) -y;
		dgdy = @(y) 2*y';

		y_true = NaN;
		z_true = NaN;
		u_true = NaN;

	elseif strcmpi(problem_str, 'alternate pendulum')
		% pendulum-spring model from eqs. (16) of Jay 2009

		g0 = 9.81; % gravitational constant
		L = 2; % length of rod
		m = 5; % mass
		k0 = 3000; % spring stiffness
		c0 = 100; % damping coefficient

		M = [m 0 0; 0 m 0; 0 0 m*L^2/3]; % mass matrix
		Minv = [1/m 0 0; 0 1/m 0; 0 0 3/(m*L^2)]; % inverse of mass matrix

		a0 = @(y,z) [0; -m*g0; -c0*z(3)-k0*(y(3)-3*pi/2)];
		b0 = @(y) [-1 0; 0 -1; -L*sin(y(3)) L*cos(y(3))];

		f = @(y,z) z;
		k = @(y,z,u) Minv*(a0(y,z)+b0(y)*u);
		g = @(y) [y(1)-L*cos(y(3)); y(2)-L*sin(y(3))];

		da0dy = [0 0 0; 0 0 0; 0 0 -k0];
		da0dz = [0 0 0; 0 0 0; 0 0 -c0];
		db0udy = @(y,u) [0 0 0; 0 0 0; 0 0 -L*u(2)*sin(y(3))-L*u(1)*cos(y(3))];

		dgdy = @(y) [1 0 L*sin(y(3)); 0 1 -L*cos(y(3))];
		dfdy = @(y,z) zeros(3);
		dfdz = @(y,z) eye(3);

		dkdy = @(y,z,u) Minv*(da0dy + db0udy(y,u));
		dkdz = @(y,z,u) Minv*da0dz;
		dkdu = @(y,z,u) Minv*b0(y);

		% the function g_{yy}(f,f)(y,z)
		gyyff = @(y,z) L*z(3)^2*[cos(y(3)); sin(y(3))];

		y30 = 3*pi/2;
		z30 = 10;

		y10 = L*cos(y30);
		y20 = L*sin(y30);

		z10 = -L*sin(y30)*z30;
		z20 = L*cos(y30)*z30;

		y0 = [y10; y20; y30];
		z0 = [z10; z20; z30];

		% Find the consistent value u0 through the find_unp function.
		% This needs u_params.  Fields not used can be set to nan.
		u_params.problem_str = problem_str;
		u_params.b = nan;
		u_params.c = nan;
		u_params.A_invertible = nan;
		u_params.Ahat_invertible = nan;
		u_params.W = nan;
		u_params.What = nan;
		u_params.nu = 2;
		u_params.g0 = g0;
		u_params.L = L;
		u_params.m = m;
		u_params.k0 = k0;
		u_params.c0 = c0;
		u_params.a0 = a0;
		u_params.b0 = b0;
		u_params.Minv = Minv;
		u_params.dgdy = dgdy;

		u0 = find_unp(y0,z0,nan,nan,u_params);

		y_true = NaN;
		z_true = NaN;
		u_true = NaN;

	elseif strcmpi(problem_str, 'nonlinear')
		% This is example 7.2 of Jay 1995

		y0 = [1; 1];
		z0 = [1; 1];
		u0 = 1;

		f1 = @(y,z) 2*y(1)*y(2)*z(1)*z(2);
		f2 = @(y,z) -y(1)*y(2)*z(2)^2;
		f = @(y,z) [f1(y,z); f2(y,z)];
		%f = @(y,z) [2*y(1)*y(2)*z(1)*z(2); -y(1)*y(2)*z(2)^2];
		k = @(y,z,u) [(y(1)*y(2)+z(1)*z(2))*u; -y(1)*y(2)^2*z(2)^3*u^2];
		g = @(y) y(1)*y(2)^2-1;

		dfdy = @(y,z) [2*y(2)*z(1)*z(2) 2*y(1)*z(1)*z(2);
						-y(2)*z(2)^2 -y(1)*z(2)^2];

		dfdz = @(y,z) [2*y(1)*y(2)*z(2) 2*y(1)*y(2)*z(1);
						0 -2*y(1)*y(2)*z(2)];

		dkdy = @(y,z,u) [y(2)*u y(1)*u;
							-y(2)^2*z(2)^3*u^2 -2*y(1)*y(2)*z(2)^3*u^2];

		dkdz = @(y,z,u) [z(2)*u z(1)*u;
							0 -3*y(1)*y(2)^2*z(2)^2*u^2];

		dkdu = @(y,z,u) [y(1)*y(2)+z(1)*z(2); -2*y(1)*y(2)^2*z(2)^3*u];

		dgdy = @(y) [y(2)^2 2*y(1)*y(2)];

		d2gdy1dy2 = @(y) 2*y(2);
		d2gdy22 = @(y) 2*y(1);

		gyyff = @(y,z) 2*d2gdy1dy2(y)*f1(y,z)*f2(y,z) + d2gdy22(y)*f2(y,z)^2;

		C = @(y,z,u) gyyff(y,z) + dgdy(y)*dfdy(y,z)*f(y,z) + dgdy(y)*dfdz(y,z)*k(y,z,u);
		dCdu = @(y,z,u) dgdy(y)*dfdz(y,z)*dkdu(y,z,u);

		y_true = @(t) [exp(2*t); exp(-t)]; 
		z_true = @(t) [exp(2*t); exp(-t)]; 
		u_true = @(t) exp(t);

	elseif strcmpi(problem_str, 'second nonlinear')
		% nonlinear problem adapted from eq. (15) of Jay 2009

		y0 = [1; 1; 0];
		z0 = [1; -2];
		u0 = 1;

		% non-constant mass matrix
		M = @(y) [y(1) y(2)-exp(-2*y(3)); sin(y(1)-exp(y(3))) y(1)*y(2)];
		%Minv = @(y) inv(M(y));

		%detM = @(y) y(1)^2*y(2) + (y(2)-exp(-2*y(3)))*sin(exp(y(3))-y(1));
		%Minv = @(y) (1/detM(y)) * [y(1)*y(2) -y(2)+exp(-2*y(3)); sin(exp(y(3))-y(1)) y(1)];

		% exact solution for M
		%M = @(y) [exp(y(3)) 0; 0 exp(-y(3))];
		%Minv = @(y) [exp(-y(3)) 0; 0 exp(y(3))];
		%M = @(y) eye(2);
		%Minv = @(y) eye(2);

		psi1 = @(t) exp(t);

		f = @(y,z) [z; 1];
		k = @(y,z,u) M(y) \ [exp(y(3))*(y(1)*z(2)+2*y(2)*z(1))+exp(2*y(3))*y(1)*u-y(1)*z(2)*psi1(y(3))-2;
			exp(-y(3))*((1/2)*y(2)*z(2)-2*y(1)*z(1)*y(2)*z(2)+y(2)*u^2)-y(1)*y(2)*z(1)*psi1(y(3))^3+exp(3*y(3))];
		g = @(y) y(1)^2*y(2)-1;

		%k = @(y,z,u) Minv(y) * [exp(y(3))*(y(1)*z(2)+2*y(2)*z(1))+exp(2*y(3))*y(1)*u-y(1)*z(2)*psi1(y(3))-2;
		%	exp(-y(3))*((1/2)*y(2)*z(2)-2*y(1)*z(1)*y(2)*z(2)+y(2)*u^2)-y(1)*y(2)*z(1)*psi1(y(3))^3+exp(3*y(3))];

		dfdy = @(y,z) zeros(3);
		dfdz = @(y,z) [eye(2); zeros(1,2)];

		dk1dy1 = @(y,z,u) exp(-5*y(3))*((-1+exp(2*y(3))*y(2))*(-4*exp(3*y(3))*y(1)*y(2)+2*exp(5*y(3))*u*y(1)^2*y(2)+2*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+exp(4*y(3))*(2+2*y(1)*y(2)*(-1+2*y(2))*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))*cos(exp(y(3))-y(1))+2*y(2)*(-exp(2*y(3))*y(1)*(-2*exp(3*y(3))*y(1)*y(2)+exp(6*y(3))*y(2)*(-2+y(1)*y(2)*z(1))+exp(4*y(3))*(2+y(1)*y(2)*(-1+2*y(2))*z(1))+y(2)*(2*u^2+z(2)-2*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-2*y(1)*z(1)*z(2)))+(-1+exp(2*y(3))*y(2))*(-2*exp(3*y(3))+2*exp(5*y(3))*u*y(1)+exp(6*y(3))*y(2)*z(1)+exp(4*y(3))*(-1+2*y(2))*z(1)-2*z(1)*z(2)+2*exp(2*y(3))*y(2)*z(1)*z(2))*sin(exp(y(3))-y(1))));

		dk1dy2 = @(y,z,u) exp(-5*y(3))*(exp(4*y(3))*y(1)^2*(2*exp(4*y(3))*y(1)*y(2)^2*z(1)+exp(2*y(3))*(-2+4*y(1)*y(2)^2*z(1))-y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))-(-4*exp(3*y(3))*y(1)+2*exp(5*y(3))*u*y(1)^2+2*u^2*(-1+exp(2*y(3))*y(2))^2-4*exp(6*y(3))*y(1)*(-1+y(2))*y(2)*z(1)-2*exp(8*y(3))*y(1)*y(2)^2*z(1)+z(2)-4*y(1)*z(1)*z(2)+2*exp(2*y(3))*y(2)*(-1+4*y(1)*z(1))*z(2)+exp(4*y(3))*(y(2)^2*z(2)-2*y(1)*z(1)*(1-4*y(2)+2*y(2)^2*z(2))))*sin(exp(y(3))-y(1)));

		dk1dy3 = @(y,z,u) exp(-5*y(3))*(exp(2*y(3))*y(1)^2*y(2)*(4*exp(5*y(3))*u*y(1)^2*y(2)+6*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+exp(4*y(3))*(2+2*y(1)*y(2)*(-1+2*y(2))*z(1))-3*y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))+exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))-exp(y(3))*(-1+exp(2*y(3))*y(2))*(-4*exp(3*y(3))*y(1)*y(2)+2*exp(5*y(3))*u*y(1)^2*y(2)+2*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+exp(4*y(3))*(2+2*y(1)*y(2)*(-1+2*y(2))*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))*cos(exp(y(3))-y(1))+(8*exp(3*y(3))*y(1)*y(2)-8*exp(5*y(3))*u*y(1)^2*y(2)+4*exp(7*y(3))*u*y(1)^2*y(2)^2+6*exp(8*y(3))*y(2)^2*(-1+y(1)*y(2)*z(1))+4*exp(6*y(3))*y(2)*(3+y(1)*(-3+y(2))*y(2)*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-2*exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2))+exp(4*y(3))*(-6+2*u^2*y(2)^3+y(2)^3*z(2)-2*y(1)*y(2)*z(1)*(-3+6*y(2)+2*y(2)^2*z(2))))*sin(exp(y(3))-y(1)));

		dk2dy1 = @(y,z,u) exp(-3*y(3))*(-exp(2*y(3))*y(1)^2*y(2)*(2*exp(4*y(3))+y(2)*(2*u^2+z(2)))-y(1)*(-4*exp(3*y(3))*y(1)*y(2)+2*exp(5*y(3))*u*y(1)^2*y(2)+2*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+exp(4*y(3))*(2+2*y(1)*y(2)*(-1+2*y(2))*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))*cos(exp(y(3))-y(1))+(8*exp(3*y(3))*y(1)*y(2)-2*exp(5*y(3))*u*y(1)^2*y(2)-2*exp(6*y(3))*y(2)*(-1+2*y(1)*y(2)*z(1))+exp(4*y(3))*(-2+4*y(1)*(1-2*y(2))*y(2)*z(1))-y(2)*(2*u^2+z(2)-8*y(1)*z(1)*z(2))+exp(2*y(3))*y(2)^2*(2*u^2+z(2)-8*y(1)*z(1)*z(2)))*sin(exp(y(3))-y(1))+2*exp(3*y(3))*u*(-1+exp(2*y(3))*y(2))*sin(exp(y(3))-y(1))^2);

		dk2dy2 = @(y,z,u) exp(-3*y(3))*(-2*exp(6*y(3))*y(1)^3-y(1)*(2*exp(6*y(3))+2*u^2-4*exp(3*y(3))*y(1)+2*exp(5*y(3))*u*y(1)^2-2*exp(4*y(3))*y(1)*z(1)+z(2)-4*y(1)*z(1)*z(2))*sin(exp(y(3))-y(1))-2*exp(2*y(3))*(-2*exp(y(3))+exp(3*y(3))*u*y(1)+2*z(1))*sin(exp(y(3))-y(1))^2);

		dk2dy3 = @(y,z,u) exp(-3*y(3))*(-exp(2*y(3))*y(1)^3*y(2)*(6*exp(4*y(3))*(-1+y(1)*y(2)*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))+exp(y(3))*y(1)*(-4*exp(3*y(3))*y(1)*y(2)+2*exp(5*y(3))*u*y(1)^2*y(2)+2*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+exp(4*y(3))*(2+2*y(1)*y(2)*(-1+2*y(2))*z(1))+y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))*cos(exp(y(3))-y(1))+y(1)*(4*exp(5*y(3))*u*y(1)^2*y(2)-6*exp(6*y(3))*y(2)*(-1+y(1)*y(2)*z(1))+2*exp(4*y(3))*(-5+y(1)*y(2)*(5+2*y(2))*z(1))-y(2)*(2*u^2+z(2)-4*y(1)*z(1)*z(2))-exp(2*y(3))*y(2)^2*(2*u^2+z(2)-4*y(1)*z(1)*z(2)))*sin(exp(y(3))-y(1))+4*exp(y(3))*(2-2*exp(2*y(3))*u*y(1)+exp(4*y(3))*u*y(1)*y(2)-3*exp(y(3))*y(2)*z(1)+exp(3*y(3))*y(2)^2*z(1))*sin(exp(y(3))-y(1))^2);

		dkdy = @(y,z,u) (1/(2*(y(1)^2*y(2)+(-exp(-2*y(3))+y(2))*sin(exp(y(3))-y(1)))^2)) * [dk1dy1(y,z,u) dk1dy2(y,z,u) dk1dy3(y,z,u); dk2dy1(y,z,u) dk2dy2(y,z,u) dk2dy3(y,z,u)];

		dk1dz1 = @(y,z,u) 2*y(1)*y(2)*(exp(6*y(3))*y(2)+exp(4*y(3))*(-1+2*y(2))-2*z(2)+2*exp(2*y(3))*y(2)*z(2));
		dk1dz2 = @(y,z,u) y(2)*(-1+exp(2*y(3))*y(2))*(-1+4*y(1)*z(1));
		dk2dz1 = @(y,z,u) -2*exp(2*y(3))*y(1)^2*y(2)*(exp(4*y(3))+2*z(2))+4*exp(4*y(3))*y(2)*sin(exp(y(3))-y(1));
		dk2dz2 = @(y,z,u) -exp(2*y(3))*y(1)*y(2)*(-1+4*y(1)*z(1));

		dkdz = @(y,z,u) (1/(2*(exp(3*y(3))*y(1)^2*y(2) + exp(y(3))*(-1 + exp(2*y(3))*y(2))*sin(exp(y(3)) - y(1))))) * [dk1dz1(y,z,u) dk1dz2(y,z,u); dk2dz1(y,z,u) dk2dz2(y,z,u)];

		dkdu = @(y,z,u) M(y) \ [y(1)*exp(2*y(3)); 2*y(2)*exp(-y(3))*u];
		%dkdu = @(y,z,u) Minv(y)*[y(1)*exp(2*y(3)); 2*y(2)*exp(-y(3))*u];

		dgdy = @(y) [2*y(1)*y(2) y(1)^2 0];

		gyyff = @(y,z) 2*y(2)*z(1)^2 + 4*y(1)*z(1)*z(2);

		y_true = @(t) [exp(t); exp(-2*t); t];
		z_true = @(t) [exp(t); -2*exp(-2*t)];
		u_true = @(t) exp(-t);

		% second hidden manifold and its derivative with respect to u
		C = @(y,z,u) gyyff(y,z) + dgdy(y)*dfdz(y,z)*k(y,z,u);
		dCdu = @(y,z,u) [2*y(1)*y(2) y(1)^2]*(M(y) \ [y(1)*exp(2*y(3)); 2*y(2)*exp(-y(3))*u]);
		%dCdu = @(y,z,u) [2*y(1)*y(2) y(1)^2]*(Minv(y)*[y(1)*exp(2*y(3)); 2*y(2)*exp(-y(3))*u]);

	elseif strcmpi(problem_str, 'linear')

		y0 = [1; 1];
		z0 = [1; 1];
		u0 = 1;

		f = @(y,z) [2*y(1)*y(2)*z(1)*z(2); -y(1)*y(2)*z(2)^2];
		k = @(y,z,u) [(y(1)*y(2)+z(1)*z(2))*u; -y(1)*y(2)^2*z(2)^2*u];
		g = @(y) y(1)*y(2)^2-1;

		dfdy = @(y,z) [2*y(2)*z(1)*z(2) 2*y(1)*z(1)*z(2);
						-y(2)*z(2)^2 -y(1)*z(2)^2];

		dfdz = @(y,z) [2*y(1)*y(2)*z(2) 2*y(1)*y(2)*z(1);
						0 -2*y(1)*y(2)*z(2)];

		dkdy = @(y,z,u) [y(2)*u y(1)*u;
							-y(2)^2*z(2)^2*u -2*y(1)*y(2)*z(2)^2*u];

		dkdz = @(y,z,u) [z(2)*u z(1)*u;
							0 -2*y(1)*y(2)^2*z(2)*u];

		dkdu = @(y,z,u) [y(1)*y(2)+z(1)*z(2); -y(1)*y(2)^2*z(2)^2];

		dgdy = @(y) [y(2)^2 2*y(1)*y(2)];

		y_true = @(t) [exp(2*t); exp(-t)]; 
		z_true = @(t) [exp(2*t); exp(-t)]; 
		u_true = @(t) exp(t);

	else
		error('Error (get_values_problem): Invalid choice of problem_str.')
	end

	true_solution.y_true = y_true;
	true_solution.z_true = z_true;
	true_solution.u_true = u_true;

	result.y0 = y0;
	result.z0 = z0;
	result.u0 = u0;
	result.f = f;
	result.k = k;
	result.g = g;
	result.dfdy = dfdy;
	result.dfdz = dfdz;
	result.dkdy = dkdy;
	result.dkdz = dkdz; 
	result.dkdu = dkdu;
	result.dgdy = dgdy;
	result.true_solution = true_solution;

	% Additional fields for the pendulum, alternate pendulum, and second
	% nonlinear problem
	if strcmpi(problem_str, 'pendulum')
		pendulum_struct.g0 = g0;
		pendulum_struct.l = l;
		result.pendulum_struct = pendulum_struct;
	elseif strcmpi(problem_str, 'alternate pendulum')
		pendulum_struct.g0 = g0;
		pendulum_struct.L = L;
		pendulum_struct.m = m;
		pendulum_struct.k0 = k0;
		pendulum_struct.c0 = c0;
		pendulum_struct.a0 = a0;
		pendulum_struct.b0 = b0;
		pendulum_struct.Minv = Minv;
		result.pendulum_struct = pendulum_struct;

		result.gyyff = gyyff;
	elseif strcmpi(problem_str,'nonlinear') || strcmpi(problem_str, 'second nonlinear')
		result.gyyff = gyyff;
		result.C = C;
		result.dCdu = dCdu;
		if strcmpi(problem_str, 'second nonlinear')
			result.M = M;
		end
	end
end

