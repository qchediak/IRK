% init file for solve_ODE
% initial conditions, parameters, etc

% The portion of code which requires the user to specify options is
% enclosed by lines consisting of the equality symbol (=) below.

% Choose the iterative scheme, problem, method, number of stages, good
% or trivial starting approximations, and plot option.
% Uncomment the desired option.

% ================================================================
% ----------------
% Iterative scheme
% ----------------

% Possible iterative schemes for the nonlinear system 0=F(x):
%
% Newton: The standard Newton method
% Simplified Newton: Uses DF(x0) as the approximate Jacobian.

% For methods which require one nonlinear system
iterative_scheme = 'newton';
%iterative_scheme = 'simplified newton';

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

% Note that these methods will all work if
% starting_approx_str='trivial', but only some will work if
% starting_approx_str='good'.

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
%s=3; % Radau IIA, Gauss, Lobatto IIIC, Lobatto III

% -----------------------
% Starting Approximations
% -----------------------

% Either choose 'trivial' or 'good'.
%starting_approx_str = 'trivial';
starting_approx_str = 'good';

% ------------
% Plot Option
% ------------

% Plot options:
%
% pendulum: Shows how well the pendulum stays on the manifold 0=y^T y-l^2 
%	and y^T v = 0 (where v=y').  Also plots position coordates and
% 	energy as a function of time.  Can only be used if
% 	problem_str='pendulum'.
% two body: plots the position coordinates of the numerical solution
% 	over one period of the true solution trajectory.
% component: Plot each of y_j vs. t separately.
% empty string: No plots.

plot_str = 'pendulum'; 
%plot_str = 'two body';
%plot_str = 'component'; 
%plot_str = ''; 

% Plot iterations: Boolean value determining if plot(s) of the iterations
% should be made.
plot_iterations = true;

% ------------
% Movie Option
% ------------

% If problem_str='pendulum' and this option is set to true, this will
% make a movie of the position coordinates of the pendulum.  Setting
% movie = true requires that plot_str=''.

%movie = true; 
movie = false; 

% -----------------------------
% Other User-Defined Parameters
% -----------------------------

% h: Stepsize (constant).
% tol: Tolerance used by the iterative solver.
% N: Number of steps.
h = 0.1; 
tol = 1e-11;
N = 50;
% ================================================================

% Get only the two body struct, if applicable, from get_values_problem
if strcmpi(problem_str, 'two body')
	problem_result = get_values_problem_ODE(problem_str);
	two_body_struct = problem_result.two_body_struct;
end

% -------------
% Miscellaneous
% -------------

% Check that either plots or a movie is produced, not both
if movie && strlength(plot_str) > 0
	error('Error (init): movie=true requires that plot_str is empty')
end

% check that the plot and movie options make sense for the problem
if strcmpi(plot_str, 'pendulum') && ~strcmpi(problem_str, 'pendulum') && ~strcmpi(problem_str, 'pendulum version 1') && ~strcmpi(problem_str, 'pendulum version 2')
	error('Error (init): The plotting option pendulum can only be chosen for the pendulum problem.')
elseif strcmpi(plot_str, 'two body') && ~strcmpi(problem_str, 'two body')
	error('Error (init): The plotting option two body can only be chosen for the two body problem.')
end

if movie && ~strcmpi(problem_str, 'pendulum') && ~strcmpi(problem_str, 'two body')
	error('Error (init): movie=true requires that either problem_str=pendulum or problem_str=two body')
end

% parameters for make_plots and make_movie
% For the pendulum problem, the field E (for energy) will also be added
% to this array after the function call find_trajectory.
plots_params.h = h;
plots_params.N = N;
plots_params.s = s;
plots_params.problem_str = problem_str;
plots_params.method_str = method_str;
plots_params.plot_str = plot_str;
plots_params.iterative_scheme = iterative_scheme;
plots_params.plot_iterations = plot_iterations;

if strcmpi(problem_str, 'two body')
	plots_params.two_body_struct = two_body_struct;
end

