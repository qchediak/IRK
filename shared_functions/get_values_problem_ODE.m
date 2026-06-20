function result = get_values_problem_ODE(problem_str)
	% function result = get_values_problem_ODE(problem_str)
	%
	% Returns values based on problem_str.  For use when solving ODEs.
	%
	% Inputs: 
	%	problem_str: String describing the ODE to be solved.
	%
	% Outputs:
	%	result: Structured array with the following fields:
	%		y0: Initial value
	%		f: The RHS of the ODE y'=f(y)
	%		dfdy: The Jacobian of f.
	%
	%	For the two body problem, result also includes the field
	%	two_body_struct, which itself has the fields a, alpha0, and T.

	% Note that the alpha0 of the pendulum and the alpha0 of the two
	% body problem are two different quantities.

	% To try a custom problem, define y0, f, and dfdy in a separate branch
	% below.  If the true solution is known, it could also be defined here.

	if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1')

		% parameters related to the pendulum
		g0 = 9.8; % Graviational constant.
		alpha0 = [0;g0];
		l = 1; % length of rod

		% initial values.  
		y0 = [3/5; -4/5; -4/2; -3/2]; 

		M1 = [1 0; 0 -1];
		M2 = [0 1; 1 0];
		M3 = [0 1; 0 0];

		f1 = @(w,z) z;
		f2 = @(w,z) (-1/dot(w,w))*(dot(z,z)-dot(w,alpha0))*w-alpha0;

		f = @(y) [f1(y(1:2),y(3:4)); f2(y(1:2),y(3:4))];

		df1dw = @(w,z) zeros(2);
		df1dz = @(w,z) eye(2);

		% Note that the block elements leave off a factor of 1/(w^Tw)^2
		df21dw1 = @(w,z) (w'*M1*w)*(dot(z,z) - dot(w,alpha0));
		df21dw2 = @(w,z) (w'*M3*alpha0)*(w'*M1*w) + (w'*M2*w)*dot(z,z);
		df22dw1 = @(w,z) (w'*M2*w)*(dot(z,z) - dot(w,alpha0));
		df22dw2 = @(w,z) (w'*M3*alpha0)*(w'*M2*w) + (w'*M1*w)*dot(z,z);

		df2dw = @(w,z) (1 / dot(w,w)^2) * [df21dw1(w,z) df21dw2(w,z); df22dw1(w,z) df22dw2(w,z)];
		df2dz = @(w,z) (-2 / dot(w,w)) * kron(w,z');

		dfdy = @(y) [df1dw(y(1:2),y(3:4)) df1dz(y(1:2),y(3:4)); 
						df2dw(y(1:2),y(3:4)) df2dz(y(1:2),y(3:4))]; 

		y_true = NaN;

	elseif strcmpi(problem_str, 'pendulum version 2')

		% parameters related to the pendulum
		g0 = 9.8; % Graviational constant.
		alpha0 = [0;g0];
		l = 1; % length of rod

		% initial values.  
		y01 = [3/5; -4/5]; 
		y02 = [-4/2; -3/2];
		y03 = (1/l^2) * (dot(y02,y02) - dot(y01,alpha0)); 

		y0 = [y01; y02; y03];

		f1 = @(y) y(3:4);
		f2 = @(y) -y(5)*y(1:2) - alpha0;
		f3 = @(y) (1/(2*dot(y(1:2),y(1:2)))) * (-8*y(5)*dot(y(1:2),y(3:4))-6*dot(y(3:4),alpha0));

		f = @(y) [f1(y); f2(y); f3(y)];

		J11 = @(y) zeros(2);
		J12 = @(y) eye(2);
		J13 = @(y) zeros(2,1);

		J21 = @(y) -y(5)*eye(2);
		J22 = @(y) zeros(2);
		J23 = @(y) -y(1:2);

		J31 = @(y) (1 / dot(y(1:2),y(1:2))^2) * ( (8*y(5)*dot(y(1:2),y(3:4)) + 6*dot(y(3:4),alpha0)) * y(1:2)' - (4*dot(y(1:2),y(1:2))*y(5)) * y(3:4)' );
		J32 = @(y) (-1 / dot(y(1:2),y(1:2))) * (4*y(5)*y(1:2)' + 3*alpha0');
		J33 = @(y) (-1 / dot(y(1:2),y(1:2))) * 4 * dot(y(1:2),y(3:4));
		
		dfdy = @(y) [J11(y) J12(y) J13(y); J21(y) J22(y) J23(y); J31(y) J32(y) J33(y)];

	elseif strcmpi(problem_str, 'two body')

		e = 0.6; % eccentricity of the orbit

		% starting values 
		% from Geometric Numerical Integration by Hairer et al, page 12
		q0 = [1-e; 0];
		p0 = [0; sqrt((1+e)/(1-e))];

		y0 = [q0; p0];

		% Semi-major axis of elliptical orbit. 
		% See, for example, "Classical Dynamics of Particles and Systems",
		% Fifth Edition, by Thornton and Marion, eq. 8.42, p. 301.
		% (The constant k=GMm in that equation is equal to 1 for our problem.)
		E0 = (1/2)*norm(p0)^2 - 1/norm(q0); % initial energy
		a = 1 / (2*abs(E0));  % semi-major axis 

		% A constant used in the plot of the true solution trajectory.
		% Ibid, eq. 8.42, p. 301
		alpha0 = a*(1-e^2);

		% Orbital period in terms of a.
		% See, for example, "Classical Dynamics of Particles and Systems",
		% Fifth Edition, by Thornton and Marion, eq. 8.49, p. 303.
		% (The expression GM in that equation is equal to 1 for our problem.)
		T = 2*pi*a^(3/2); 
 
		f1 = @(y) y(3:4);
		f2 = @(y) -(dot(y(1:2),y(1:2))^(-3/2))*y(1:2);

		f = @(y) [f1(y); f2(y)];

		J11 = @(y) zeros(2);
		J12 = @(y) eye(2);
		J21 = @(y) 3*(dot(y(1:2),y(1:2))^(-5/2))*(y(1:2)*y(1:2)') - (dot(y(1:2),y(1:2))^(-3/2))*eye(2);
		J22 = @(y) zeros(2);

		dfdy = @(y) [J11(y) J12(y); J21(y) J22(y)];

	else

		error('Error (get_values_problem): Invalid choice of problem_str.')

	end

	result.y0 = y0;
	result.f = f;
	result.dfdy = dfdy;

	% Additional fields for specific problems
	if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1') || strcmpi(problem_str, 'pendulum version 2')
		pendulum_struct.g0 = g0;
		pendulum_struct.l = l;

		result.pendulum_struct = pendulum_struct;
	elseif strcmpi(problem_str, 'two body')
		two_body_struct.a = a;
		two_body_struct.alpha0 = alpha0;
		two_body_struct.T = T;
		two_body_struct.e = e;

		result.two_body_struct = two_body_struct;
	end
end

