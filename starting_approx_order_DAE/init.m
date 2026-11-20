% init.m
% initial conditions, parameters, etc 

% The portion of code which requires the user to specify options is
% enclosed by lines consisting of the equality symbol (=) below.

% Choose the iterative scheme, problem, method, number of stages
% and whether the starting approximations should be good
% or trivial. Uncomment the desired option.

% ================================================================
% ----------------
% Iterative scheme
% ----------------

% Possible iterative schemes for the nonlinear system 0=F(x):
%
% Newton: The standard Newton method
% Simplified Newton: Uses DF(x0) as the approximate Jacobian.
% Simplified Newton 0: Uses constraint equation without scaling.
% Simplified Newton 1: Same as simplified newton.  The "1" is used in
% 	some code to contrast with the "2" option below.
% Simplified Newton 2: Approximates DF(x0) by another matrix J, and uses
% 	J as the approximate Jacobian. 
% Simplified Newton 3: Similar to simplified newton 2, but approximates
% 	more entries of J as zero.
% kelley: Uses a more sophisticated root finder (nsold) by C.T. Kelley.
% 	Available here: https://github.com/ctkelley/ctkelley.github.io/blob/b0aa502875f9409a16254d81de37a8d9a7ce0cf4/newton/SOLVERS/nsold.m


% For methods which require one nonlinear system
%iterative_scheme = 'newton';
%iterative_scheme = 'simplified newton 1';
%iterative_scheme = 'simplified newton 2'; 
%iterative_scheme = 'simplified newton 3'; 
%iterative_scheme = 'kelley';

% For methods which require two nonlinear systems
%iterative_scheme1 = 'newton';
%iterative_scheme1 = 'simplified newton 0'; % Lobatto IIIA-IIIB only
iterative_scheme1 = 'simplified newton';
%iterative_scheme1 = 'simplified newton 2';
%iterative_scheme1 = 'simplified newton 3';
%iterative_scheme1 = 'kelley';
%iterative_scheme2 = 'newton';
iterative_scheme2 = 'simplified newton';
%iterative_scheme2 = 'kelley';

% -------
% Problem
% -------

% pendulum: The plane pendulum problem.  k is linear in u.
% linear: An index 3 DAE in which k is linear in u.
% nonlinear: An example where k is nonlinear in u.  Requires that
% 	plot_str be set to component.

problem_str = 'Pendulum';
%problem_str = 'alternate pendulum';
%problem_str = 'Linear';
%problem_str = 'Nonlinear';

% -------
% Method
% -------

% Note that there is a separate method_str for Lobatto IIIC and
% projected Lobatto IIIC.

%method_str = 'Gauss';
%method_str = 'Radau IIA';
%method_str = 'Lobatto IIIC';
method_str = 'Lobatto IIIA-IIIB';
%method_str = 'Projected Lobatto IIIC';
%method_str = 'Projected Radau IIA';

% ----------------
% Number of stages
% ----------------

%s=2; % Radau IIA, Lobatto IIIA-IIIB
s=3; % Gauss, Lobatto IIIC, projected Lobatto IIIC, projected Radau IIA, Lobatto IIIA-IIIB

% -----------------------
% Starting Approximations
% -----------------------

% Choose one of the following:
%	'trivial' (the vector (xn,xn,...,xn)), 
%	'good' (the starting approximations based on my thesis)
%	'interpolation' (starting approximations found through the
%		interpolating polynomial).  This currently requires that s=3 and
%		tnm=0; i.e. the interpolating polynomial interpolates over the
%		interval [tnm,tn]=[0,h] and is extrapolated to the interval
%		[tn,tnp]=[h,(1+r)h].
%starting_approx_str = 'trivial';
%starting_approx_str = 'good';
starting_approx_str = 'interpolation';

% ------------
% Plot Option
% ------------

% Options for plotting errors:
%
% 1: Useful for testing.  Allows more flexibility in what plots to make,
% 	such as 'component' or 'composite' (see below), and separately plots
% 	the errors of znp0, Unps0, or Unpsp0 if applicable.
%
% 2: The option used to make some of the graphs in my thesis.  
%	Requires the method be projected Lobatto IIIC with s=3.
%	Plots the y-components separately, combines znp0 with Z0 and Unpsp0
%	with U0.
error_plotting_option = 1;

% Options for plotting errors: 
% Either 'composite' (i.e. error in all of X_{n+1}^{(0)})
% or 'component' (i.e. error in each X_{n+1,i}^{(0)} for i=1,...,s).
% Currently the 'component' option is only available for s=3.
% These only change the plots if error_plotting_option=1.
error_plotting_Y = 'component';
error_plotting_Z = 'component';
error_plotting_U = 'component';

% -----------------------------
% Other User-Defined Parameters
% -----------------------------

% Option for find_error_SA.
% The choices are 'BF', 'FF', and 'BFF'.  
% See find_error_SA documentation for descriptions of these.
% 'BF' is the recommended option for most methods.
% Use the option 'FF' for Radau IIA, Lobatto IIIC, and 
% projected Radau IIA.
find_error_SA_option = 'FF';

% tol: Tolerance used by the iterative solver.
% r: r=h2/h1, where h1 and h2 are the stepsizes of the first and second
% 	step, respectively.
% stepsizes_len: Number of stepsizes to test.  
% logspace_lower_lim: Lower limit of stepsizes in logspace.
% logpsace_upper_lim: Upper limit of stepsizes in logspace.
tol = 1e-10;
%r = 1+sqrt(15)/3; % special value for Gauss
r = 1/2;
stepsizes_len = 4;
logspace_lower_lim = -3.8; 
logspace_upper_lim = -3.7; 
% ================================================================

stepsizes = logspace(logspace_lower_lim, logspace_upper_lim, stepsizes_len);
stepsizes = flip(stepsizes);

% check that error_plotting_option makes sense
if error_plotting_option==2 && (~strcmpi(method_str,'projected Lobatto IIIC') || s~=3)
	error('Error (init): error_plotting_option=2 requires the method be 3-stage projected Lobatto IIIC')
end

% ----------------------
% Configure the method
% ----------------------

% Get values based on method_str and s
method_result = get_values_method(method_str, s);

% unpack method_result
A_invertible = method_result.A_invertible;
Ahat_invertible = method_result.Ahat_invertible;
Astar_invertible = method_result.Astar_invertible;
Ahatstar_invertible = method_result.Ahatstar_invertible;
W = method_result.W;
What = method_result.What;
b = method_result.b;
c = method_result.c;
A = method_result.A;
Ahat = method_result.Ahat;
bstar = method_result.bstar;
cstar = method_result.cstar;
Astar = method_result.Astar;
Ahatstar = method_result.Ahatstar;
Wstar = method_result.Wstar;
Whatstar = method_result.Whatstar;
number_nonlinear_systems = method_result.number_nonlinear_systems;
reverse_method_str = method_result.reverse_method_str;

% ----------------------
% Configure the problem
% ----------------------

problem_result = get_values_problem(problem_str);

% unpack problem_result
y0 = problem_result.y0;
z0 = problem_result.z0;
u0 = problem_result.u0;
f = problem_result.f;
k = problem_result.k;
g = problem_result.g;
dfdy = problem_result.dfdy;
dfdz = problem_result.dfdz;
dkdy = problem_result.dkdy;
dkdz = problem_result.dkdz; 
dkdu = problem_result.dkdu;
dgdy = problem_result.dgdy;
true_solution = problem_result.true_solution;

y_true = true_solution.y_true;
z_true = true_solution.z_true;
u_true = true_solution.u_true;

% Additional values for the pendulum
if strcmpi(problem_str, 'pendulum')
	g0 = problem_result.pendulum_struct.g0;
	l = problem_result.pendulum_struct.l;
	alpha0 = [0;g0];
elseif strcmpi(problem_str, 'alternate pendulum')
	g0 = problem_result.pendulum_struct.g0;
	L = problem_result.pendulum_struct.L;
	m = problem_result.pendulum_struct.m;
	k0 = problem_result.pendulum_struct.k0;
	c0 = problem_result.pendulum_struct.c0;
	a0 = problem_result.pendulum_struct.a0;
	b0 = problem_result.pendulum_struct.b0;
	Minv = problem_result.pendulum_struct.Minv;
elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
	gyyff = problem_result.gyyff;
	C = problem_result.C;
	dCdu = problem_result.dCdu;
	if strcmpi(problem_str, 'second nonlinear')
		M = problem_result.M;
	end
end

% Note that this assumes y0, z0, u0 are column vectors
ny = size(y0,1); % dimensions of y
nz = size(z0,1); % dimensions of z
nu = size(u0,1); % dimensions of u

% --------------------
% Iterative Scheme(s)
% --------------------

% Check that the choice of iterative_scheme1 and iterative_scheme2 makes sense
if strcmpi(iterative_scheme1,'simplified newton 0') && ~strcmpi(method_str,'Lobatto IIIA-IIIB')
	error('Error (init): The iterative scheme simplified newton 0 can only be used with the Lobatto IIIA-IIIB method')
elseif strcmpi(iterative_scheme2,'simplified newton 2') 
	error('Error (init): The iterative scheme simplified newton 2 cannot be used with the second nonlinear system.')
end

% Define iterative_scheme_struct
if number_nonlinear_systems==1
	iterative_scheme_struct.iterative_scheme = iterative_scheme;
elseif number_nonlinear_systems==2
	iterative_scheme_struct.iterative_scheme1 = iterative_scheme1;
	iterative_scheme_struct.iterative_scheme2 = iterative_scheme2;
else
	error('Error (init): Invalid value of number_nonlinear_systems')
end

% --------------
% Miscellaneous
% --------------

% Check that find_error_SA_option makes sense with the given method.
% In particular, 'BF' and 'BFF' cannot be used with projected Radau IIA
% or Lobatto IIIC, because the reverse methods (projected Radau I and
% Lobatto III, respectively) are not implemented here.
% For these methods, use 'FF'.
if (strcmpi(method_str, 'projected Radau IIA') || strcmpi(method_str, 'Radau IIA') || strcmpi(method_str, 'Lobatto IIIC')) && (strcmpi(find_error_SA_option, 'BF') || strcmpi(find_error_SA_option, 'BFF'))
	error('Error (init): Invalid combination of method_str and find_error_SA_option.  Use the option FF instead.')
end


% parameters for find_error_SA
params.starting_approx_str = starting_approx_str;
params.problem_str = problem_str;
params.method_str = method_str;
params.iterative_scheme_struct = iterative_scheme_struct;
params.r = r;
params.f = f;
params.k = k;
params.g = g;
params.dfdy = dfdy;
params.dfdz = dfdz;
params.dkdy = dkdy;
params.dkdz = dkdz;
params.dkdu = dkdu; 
params.dgdy = dgdy;
params.ny = ny;
params.nz = nz;
params.nu = nu;
params.b = b;
params.c = c;
params.A = A;
params.Ahat = Ahat;
params.bstar = bstar;
params.cstar = cstar;
params.Astar = Astar;
params.Ahatstar = Ahatstar;
params.tol = tol;
params.y0 = y0;
params.z0 = z0;
params.u0 = u0;
params.A_invertible = A_invertible;
params.Ahat_invertible = Ahat_invertible;
params.Astar_invertible = Astar_invertible;
params.Ahatstar_invertible = Ahatstar_invertible;
params.W = W;
params.What = What;
params.Wstar = Wstar;
params.Whatstar = Whatstar;
params.number_nonlinear_systems = number_nonlinear_systems;
params.reverse_method_str = reverse_method_str;

if strcmpi(problem_str, 'pendulum')
	pendulum_params.g0 = g0;
	pendulum_params.l = l;
	params.pendulum_params = pendulum_params;
elseif strcmpi(problem_str, 'alternate pendulum')
	pendulum_params.g0 = g0;
	pendulum_params.L = L;
	pendulum_params.m = m;
	pendulum_params.k0 = k0;
	pendulum_params.c0 = c0;
	pendulum_params.a0 = a0;
	pendulum_params.b0 = b0;
	pendulum_params.Minv = Minv;
	params.pendulum_params = pendulum_params;
elseif strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'second nonlinear')
	nonlinear_params.gyyff = gyyff;
	nonlinear_params.C = C;
	nonlinear_params.dCdu = dCdu;
	if strcmpi(problem_str, 'second nonlinear')
		nonlinear_params.M = M;
	end

	params.nonlinear_params = nonlinear_params;
end

% Structured array of error plotting options
error_plotting_struct.error_plotting_option = error_plotting_option;
error_plotting_struct.error_plotting_Y = error_plotting_Y;
error_plotting_struct.error_plotting_Z = error_plotting_Z;
error_plotting_struct.error_plotting_U = error_plotting_U;

% parameters for printing and plotting
print_plot_params.problem_str = problem_str;
print_plot_params.method_str = method_str;
print_plot_params.starting_approx_str = starting_approx_str;
print_plot_params.number_nonlinear_systems = number_nonlinear_systems;
print_plot_params.r = r;
print_plot_params.s = s;
print_plot_params.tol = tol;
print_plot_params.iterative_scheme_struct = iterative_scheme_struct;
print_plot_params.error_plotting_struct = error_plotting_struct;

