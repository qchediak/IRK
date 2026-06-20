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
% Simplified Newton 1: Same as simplified newton.  The "1" is used in
% 	some code to contrast with the "2" option below.

iterative_scheme = 'newton';
%iterative_scheme = 'simplified newton 1';

% -------
% Problem
% -------

% 'pendulum' or 'pendulum version 1': The plane pendulum problem as an ODE.  
%	This is the same as the partitioned ODE with an equation for y' and z',
%	with the entire RHS forming f(y).
%
% 'pendulum version 2': The plane pendulum as an ODE.  This version
%	has an equation for y', z', and u', with the entire RHS forming f(y).
%
% 'two body': The two body problem as an ODE.

%problem_str = 'pendulum';
%problem_str = 'pendulum version 2';
problem_str = 'two body';

% -------
% Method
% -------

% note that not every method below has good starting approximations set

%method_str = 'Gauss';
method_str = 'Radau IIA';
%method_str = 'Lobatto IIIA';
%method_str = 'Lobatto IIIB';
%method_str = 'Lobatto IIIC';
%method_str = 'Lobatto III';

% ----------------
% Number of stages
% ----------------

s=2; % Radau IIA, Gauss, Lobatto IIIA, Lobatto IIIB
%s=3; % Gauss, Lobatto IIIC 

% -----------------------
% Starting Approximations
% -----------------------

% Either choose 'trivial' or 'good'.
%starting_approx_str = 'trivial';
starting_approx_str = 'good';

% ------------
% Plot Option
% ------------

% Options for plotting errors: 
% Either 'composite' (i.e. error in all of X_{n+1}^{(0)})
% or 'component' (i.e. error in each X_{n+1,i}^{(0)} for i=1,...,s).
%error_plotting_Y = 'component';
error_plotting_Y = 'composite';

% -----------------------------
% Other User-Defined Parameters
% -----------------------------

% Option for find_error_SA.
% The choices are 'BF', 'FF', and 'BFF'.  
% See find_error_SA documentation for descriptions of these.
% 'BF' is the recommended option for most methods.
% Use the option 'FF' for Radau IIA and Lobatto IIIC
find_error_SA_option = 'FF';

% tol: Tolerance used by the iterative solver.
% r: r=h2/h1, where h1 and h2 are the stepsizes of the first and second
% 	step, respectively.
% stepsizes_len: Number of stepsizes to test.  
% logspace_lower_lim: Lower limit of stepsizes in logspace.
% logpsace_upper_lim: Upper limit of stepsizes in logspace.
tol = 1e-11;
r = 1/2;
%r = 1+sqrt(3); % special value for Gauss (ODE)
stepsizes_len = 10;
logspace_lower_lim = -2; 
logspace_upper_lim = -4; 
% ================================================================

stepsizes = logspace(logspace_lower_lim, logspace_upper_lim, stepsizes_len);
stepsizes = flip(stepsizes);

% ----------------------
% Configure the method
% ----------------------

% Get values based on method_str and s
method_result = get_values_method(method_str, s);

% unpack method_result
A_invertible = method_result.A_invertible;
Astar_invertible = method_result.Astar_invertible;
W = method_result.W;
b = method_result.b;
c = method_result.c;
A = method_result.A;
bstar = method_result.bstar;
cstar = method_result.cstar;
Astar = method_result.Astar;
Wstar = method_result.Wstar;
reverse_method_str = method_result.reverse_method_str;

% ----------------------
% Configure the problem
% ----------------------

problem_result = get_values_problem_ODE(problem_str);

% unpack problem_result
y0 = problem_result.y0;
f = problem_result.f;
dfdy = problem_result.dfdy;

% Additional values for the pendulum
if strcmpi(problem_str, 'pendulum')
	g0 = problem_result.pendulum_struct.g0;
	l = problem_result.pendulum_struct.l;
	alpha0 = [0;g0];
elseif strcmpi(problem_str, 'two body')
	a = problem_result.two_body_struct.a;
	alpha0 = problem_result.two_body_struct.alpha0;
	T = problem_result.two_body_struct.T;
	e = problem_result.two_body_struct.e;
end

% Note that this assumes y0 is a column vector
ny = size(y0,1); % dimension of y

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
params.problem_str = problem_str;
params.method_str = method_str;
params.iterative_scheme = iterative_scheme;
params.r = r;
params.f = f;
params.dfdy = dfdy;
params.ny = ny;
params.b = b;
params.c = c;
params.A = A;
params.bstar = bstar;
params.cstar = cstar;
params.Astar = Astar;
params.tol = tol;
params.y0 = y0;
params.A_invertible = A_invertible;
params.Astar_invertible = Astar_invertible;
params.W = W;
params.Wstar = Wstar;
params.reverse_method_str = reverse_method_str;

if strcmpi(problem_str, 'pendulum')
	pendulum_params.g0 = g0;
	pendulum_params.l = l;
	params.pendulum_params = pendulum_params;
elseif strcmpi(problem_str, 'two body')
	two_body_params.a = a;
	two_body_params.alpha0 = alpha0;
	two_body_params.T = T;
	two_body_params.e = e;
	params.two_body_params = two_body_params;
end

% parameters for printing and plotting
print_plot_params.problem_str = problem_str;
print_plot_params.method_str = method_str;
print_plot_params.starting_approx_str = starting_approx_str;
print_plot_params.r = r;
print_plot_params.s = s;
print_plot_params.tol = tol;
print_plot_params.iterative_scheme = iterative_scheme;
print_plot_params.error_plotting_Y = error_plotting_Y;


