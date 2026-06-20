% init (starting_approx_order_pODE)

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

iterative_scheme = 'newton';
%iterative_scheme = 'simplified newton';

% -------
% Problem
% -------

% pendulum: The plane pendulum problem as a partitioned ODE.
% two body: The two body problem as a partitioned ODE.
% Lotka-Volterra: The Lotka-Volterra problem

%problem_str = 'pendulum';
%problem_str = 'two body';
problem_str = 'Lotka-Volterra';


% -------
% Method
% -------

% Note that these methods will all work if
% starting_approx_str='trivial', but only some will work if
% starting_approx_str='good'.

% These are all RK methods, with the exception of Lobatto IIIA-IIIB.

%method_str = 'Gauss';
%method_str = 'Radau IIA';
%method_str = 'Lobatto IIIA';
%method_str = 'Lobatto IIIB';
%method_str = 'Lobatto IIIC';
%method_str = 'Lobatto III';
method_str = 'Lobatto IIIA-IIIB';

% ----------------
% Number of stages
% ----------------

s=2; % Radau IIA, Gauss, Lobatto IIIA, Lobatto IIIB, Lobatto IIIA-IIIB
%s=3; % Radau IIA, Gauss, Lobatto IIIC, Lobatto III, Lobatto IIIA-IIIB

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
error_plotting_Y = 'component';
error_plotting_Z = 'component';

% -----------------------------
% Other User-Defined Parameters
% -----------------------------

% Option for find_error_SA.
% The choices are 'BF', 'FF', and 'BFF'.  
% See find_error_SA documentation for descriptions of these.
% 'BF' is the recommended option for most methods.
% Use the option 'FF' for Radau IIA and Lobatto IIIC
find_error_SA_option = 'BF';

% tol: Tolerance used by the iterative solver.
% r: r=h2/h1, where h1 and h2 are the stepsizes of the first and second
% 	step, respectively.
% stepsizes_len: Number of stepsizes to test.  
% logspace_lower_lim: Lower limit of stepsizes in logspace.
% logpsace_upper_lim: Upper limit of stepsizes in logspace.
tol = 1e-15;
%r = 1+sqrt(15)/3; % special value for Gauss
%r = 1/2;
r = 1;
stepsizes_len = 10;
logspace_lower_lim = -3; 
logspace_upper_lim = -3.5; 
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
reverse_method_str = method_result.reverse_method_str;

% ----------------------
% Configure the problem
% ----------------------

problem_result = get_values_problem_pODE(problem_str);

% unpack problem_result
y0 = problem_result.y0;
z0 = problem_result.z0;
f = problem_result.f;
g = problem_result.g;
dfdy = problem_result.dfdy;
dfdz = problem_result.dfdz;
dgdy = problem_result.dgdy;
dgdz = problem_result.dgdz;

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

% Note that this assumes y0 and z0 are column vectors
ny = size(y0,1); % dimensions of y
nz = size(z0,1); % dimensions of z

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
params.g = g;
params.dfdy = dfdy;
params.dfdz = dfdz;
params.dgdy = dgdy;
params.dgdz = dgdz;
params.ny = ny;
params.nz = nz;
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
params.A_invertible = A_invertible;
params.Ahat_invertible = Ahat_invertible;
params.Astar_invertible = Astar_invertible;
params.Ahatstar_invertible = Ahatstar_invertible;
params.W = W;
params.What = What;
params.Wstar = Wstar;
params.Whatstar = Whatstar;
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

% Structured array of error plotting options
error_plotting_struct.error_plotting_Y = error_plotting_Y;
error_plotting_struct.error_plotting_Z = error_plotting_Z;

% parameters for printing and plotting
print_plot_params.problem_str = problem_str;
print_plot_params.method_str = method_str;
print_plot_params.starting_approx_str = starting_approx_str;
print_plot_params.r = r;
print_plot_params.s = s;
print_plot_params.tol = tol;
print_plot_params.iterative_scheme = iterative_scheme;
print_plot_params.error_plotting_struct = error_plotting_struct;

