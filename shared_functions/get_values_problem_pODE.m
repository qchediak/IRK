function result = get_values_problem_pODE(problem_str)
	%function result = get_values_problem_pODE(problem_str)
	%
	% Returns values based on problem_str.  For use when solving partitioned
	% ODEs.
	%
	% Inputs: 
	%	problem_str: String describing the partitioned ODE to be solved.
	%
	% Outputs:
	%	result: Structured array with the following fields:
	%		y0, z0: Initial values
	%		f, g: The RHS fuctions of the partitioned ODE y'=f(y,z), z'=g(y,z)
	%		dfdy, dfdz, dgdy, dgdz: partial derivatives of f and g.
	%
	% For the two body problem, result also includes the field
	% two_body_struct, which itself has the fields a, alpha0, and T.

	% Note that the alpha0 of the pendulum and the alpha0 of the two
	% body problem are two different quantities.

	% A custom problem could be defined below. 

	if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1')

		% parameters related to the pendulum
		g0 = 9.8; % Graviational constant.
		alpha0 = [0;g0];
		l = 1; % length of rod

		% initial values.  
		y0 = [3/5; -4/5]; 
		z0 = [-4/2; -3/2];

		M1 = [1 0; 0 -1];
		M2 = [0 1; 1 0];
		M3 = [0 1; 0 0];

		f = @(y,z) z;
		g = @(y,z) (-1/dot(y,y))*(dot(z,z)-dot(y,alpha0))*y-alpha0;
		% simplified with y^Ty=l^2 
		%g = @(y,z) (-1/l^2)*(dot(z,z)-dot(y,alpha0))*y-alpha0;

		dfdy = @(y,z) zeros(2);
		dfdz = @(y,z) eye(2);

		% Note that the block elements of dgdy leave off a factor of 1/(y^Ty)^2
		dg1dy1 = @(y,z) (y'*M1*y)*(dot(z,z) - dot(y,alpha0));
		dg1dy2 = @(y,z) (y'*M3*alpha0)*(y'*M1*y) + (y'*M2*y)*dot(z,z);
		dg2dy1 = @(y,z) (y'*M2*y)*(dot(z,z) - dot(y,alpha0));
		dg2dy2 = @(y,z) (y'*M3*alpha0)*(y'*M2*y) + (y'*M1*y)*dot(z,z);

		dgdy = @(y,z) (1 / dot(y,y)^2) * [dg1dy1(y,z) dg1dy2(y,z); dg2dy1(y,z) dg2dy2(y,z)];

		dgdz = @(y,z) (-2 / dot(y,y)) * kron(y,z');

		y_true = NaN;
		z_true = NaN;

	elseif strcmpi(problem_str, 'pendulum theta') 

		% The ODE is theta'' + k0^2*sin(theta) = 0
		% The pODE is theta'=omega, omega'=-k0^2*sin(theta)

		l = 1; % length of rod
		g0 = 9.8; % gravitational constant
		k0 = sqrt(g0/l);

		y0 = atan(-4/5);		% aka theta_0
		z0 = -2 / (l*cos(y0));	% aka omega_0

		f = @(theta,omega) omega;
		g = @(theta,omega) -k0^2*sin(theta);

		dfdy = @(theta,omega) 0; % aka df/d(theta)
		dfdz = @(theta,omega) 1; % aka df/d(omega)
		dgdy = @(theta,omega) -k0^2*cos(theta); % aka dg/d(theta)
		dgdz = @(theta,omega) 0; % aka dg/d(omega)

	elseif strcmpi(problem_str, 'double pendulum angles')

		l1=1;
		l2=1;
		g0=9.8;

		y0 = [pi/2; pi]; % aka theta0
		z0 = [0; 0]; % aka omega0

		f = @(theta,omega) omega;

		dfdy = @(theta,omega) zeros(2,2); % aka df/d(theta)

		dfdz = @(theta,omega) eye(2); % aka df/d(omega)

		g1 = @(theta,omega) (3*g0*sin(theta(1))+g0*sin(theta(1)-2*theta(2))+2*(l2*omega(2)^2+l1*omega(1)^2*cos(theta(1)-theta(2)))*sin(theta(1)-theta(2)))/(l1*(-3+cos(2*(theta(1)-theta(2)))));

		g2 = @(theta,omega) -((2*(2*l1*omega(1)^2+2*g0*cos(theta(1))+l2*omega(2)^2*cos(theta(1)-theta(2)))*sin(theta(1)-theta(2)))/(l2*(-3+cos(2*(theta(1)-theta(2))))));

		g = @(theta,omega) [g1(theta,omega); g2(theta,omega)];

		dgdt11 = @(theta,omega) -((-4*l1*omega(1)^2+15*g0*cos(theta(1))+g0*cos(3*theta(1)-4*theta(2))-3*g0*cos(theta(1)-2*theta(2))+3*g0*cos(3*theta(1)-2*theta(2))+6*l2*omega(2)^2*cos(theta(1)-theta(2))+12*l1*omega(1)^2*cos(2*(theta(1)-theta(2)))+2*l2*omega(2)^2*cos(3*(theta(1)-theta(2))))/(2*l1*(-3+cos(2*(theta(1)-theta(2))))^2));

		dgdt12 = @(theta,omega) (-2*l1*omega(1)^2-2*g0*cos(theta(1))+3*g0*cos(theta(1)-2*theta(2))+3*g0*cos(3*theta(1)-2*theta(2))+3*l2*omega(2)^2*cos(theta(1)-theta(2))+6*l1*omega(1)^2*cos(2*(theta(1)-theta(2)))+l2*omega(2)^2*cos(3*(theta(1)-theta(2))))/(l1*(-3+cos(2*(theta(1)-theta(2))))^2);

		dgdt21 = @(theta,omega) (2*(-l2*omega(2)^2+g0*cos(2*theta(1)-3*theta(2))+3*l1*omega(1)^2*cos(theta(1)-theta(2))+3*l2*omega(2)^2*cos(2*(theta(1)-theta(2)))+l1*omega(1)^2*cos(3*(theta(1)-theta(2)))+5*g0*cos(2*theta(1)-theta(2))-2*g0*cos(theta(2))))/(l2*(-3+cos(2*(theta(1)-theta(2))))^2);

		dgdt22 = @(theta,omega) -((-2*l2*omega(2)^2+g0*cos(2*theta(1)-3*theta(2))+g0*cos(4*theta(1)-3*theta(2))+6*l1*omega(1)^2*cos(theta(1)-theta(2))+6*l2*omega(2)^2*cos(2*(theta(1)-theta(2)))+2*l1*omega(1)^2*cos(3*(theta(1)-theta(2)))+3*g0*cos(2*theta(1)-theta(2))+3*g0*cos(theta(2)))/(l2*(-3+cos(2*(theta(1)-theta(2))))^2));

		dgdw11 = @(theta,omega) (2*omega(1)*sin(2*(theta(1)-theta(2))))/(-3+cos(2*(theta(1)-theta(2))));

		dgdw12 = @(theta,omega) (4*l2*omega(2)*sin(theta(1)-theta(2)))/(l1*(-3+cos(2*(theta(1)-theta(2)))));

		dgdw21 = @(theta,omega) -((8*l1*omega(1)*sin(theta(1)-theta(2)))/(l2*(-3+cos(2*(theta(1)-theta(2))))));

		dgdw22 = @(theta,omega) -((2*omega(2)*sin(2*(theta(1)-theta(2))))/(-3+cos(2*(theta(1)-theta(2)))));

		dgdy = @(theta,omega) [dgdt11(theta,omega) dgdt12(theta,omega);
			dgdt21(theta,omega) dgdt22(theta,omega)];

		dgdz = @(theta,omega) [dgdw11(theta,omega) dgdw12(theta,omega);
			dgdw21(theta,omega) dgdw22(theta,omega)];

	elseif strcmpi(problem_str, 'linear')

		% note that t=z3 to make this autonomous
		y0=[1;1];
		z0=[1;1;0];

		f = @(y,z) [2*y(1)*y(2)*z(1)*z(2);
			-y(1)*y(2)*z(2)^2];

		g = @(y,z) [(y(1)*y(2)+z(1)*z(2))*exp(z(3));
			-y(1)*y(2)^2*z(2)^2*exp(z(3));
			1];

		dfdy = @(y,z) [2*y(2)*z(1)*z(2) 2*y(1)*z(1)*z(2);
			-y(2)*z(2)^2 -y(1)*z(2)^2];

		dfdz = @(y,z) [2*y(1)*y(2)*z(2) 2*y(1)*y(2)*z(1) 0;
			-y(2)*z(2)^2 -y(1)*z(2)^2 0];

		dgdy = @(y,z) [y(2)*exp(z(3)) y(1)*exp(z(3));
			-y(2)^2*z(2)^2*exp(z(3)) -2*y(1)*y(2)*z(2)^2*exp(z(3));
			0 0];

		dgdz = @(y,z) [z(2)*exp(z(3)) z(1)*exp(z(3)) (y(1)*y(2)+z(1)*z(2))*exp(z(3));
			0 -2*y(1)*y(2)^2*z(2)*exp(z(3)) -y(1)*y(2)^2*z(2)^2*exp(z(3));
			0 0 0];

	elseif strcmpi(problem_str, 'two body')

		e = 0.6; % eccentricity of the orbit

		% starting values 
		% from Geometric Numerical Integration by Hairer et al, page 12
		y0 = [1-e; 0];
		z0 = [0; sqrt((1+e)/(1-e))];

		% Semi-major axis of elliptical orbit. 
		% See, for example, "Classical Dynamics of Particles and Systems",
		% Fifth Edition, by Thornton and Marion, eq. 8.42, p. 301.
		% (The constant k=GMm in that equation is equal to 1 for our problem.)
		E0 = (1/2)*norm(z0)^2 - 1/norm(y0); % initial energy
		a = 1 / (2*abs(E0));  % semi-major axis 

		% A constant used in the plot of the true solution trajectory.
		% Ibid, eq. 8.42, p. 301
		alpha0 = a*(1-e^2);

		% Orbital period in terms of a.
		% See, for example, "Classical Dynamics of Particles and Systems",
		% Fifth Edition, by Thornton and Marion, eq. 8.49, p. 303.
		% (The expression GM in that equation is equal to 1 for our problem.)
		T = 2*pi*a^(3/2); 
 
		f = @(y,z) z;
		g = @(y,z) -(dot(y,y)^(-3/2))*y;

		dfdy = @(y,z) zeros(2);
		dfdz = @(y,z) eye(2);

		dgdy = @(y,z) 3*(dot(y,y)^(-5/2))*(y*y') - (dot(y,y)^(-3/2))*eye(2);
		dgdz = @(y,z) zeros(2);

	elseif strcmpi(problem_str, 'lotka-volterra')

		% initial values
		y0 = 2;
		z0 = 2;

		% constants
		c0=2;
		a0=1;
		gamma0=1;
		alpha0=1;

		%f = @(y,z) y*(z-2);
		%g = @(y,z) z*(1-y);
		%dfdy = @(y,z) z-2;
		%dfdz = @(y,z) y;
		%dgdy = @(y,z) -z;
		%dgdz = @(y,z) 1-y;

		f = @(y,z) -c0*y + gamma0*y*z;
		g = @(y,z) a0*z - alpha0*y*z;

		dfdy = @(y,z) -c0 + gamma0*z;
		dfdz = @(y,z) gamma0*y;
		dgdy = @(y,z) -alpha0*z;
		dgdz = @(y,z) a0 - alpha0*y;

	else

		error('Error (get_values_problem): Invalid choice of problem_str.')

	end

	result.y0 = y0;
	result.z0 = z0;
	result.f = f;
	result.g = g;
	result.dfdy = dfdy;
	result.dfdz = dfdz;
	result.dgdy = dgdy;
	result.dgdz = dgdz;

	% Additional fields for specific problems
	if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum theta')
		pendulum_struct.g0 = g0;
		pendulum_struct.l = l;

		result.pendulum_struct = pendulum_struct;
	elseif strcmpi(problem_str, 'double pendulum angles')
		pendulum_struct.g0 = g0;
		pendulum_struct.l1 = l1;
		pendulum_struct.l2 = l2;

		result.pendulum_struct = pendulum_struct;
	elseif strcmpi(problem_str, 'two body')
		two_body_struct.a = a;
		two_body_struct.alpha0 = alpha0;
		two_body_struct.T = T;
		two_body_struct.e = e;

		result.two_body_struct = two_body_struct;
	end
end

