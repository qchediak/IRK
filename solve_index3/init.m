% init.m
% initial conditions, parameters, etc

% The portion of code which requires the user to specify options is
% enclosed by lines consisting of the equality symbol (=) below.

% Choose the iterative scheme, problem, method, number of stages, good
% or trivial starting approximations, and plot and movie option.
% Uncomment the desired option.

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
% kelley: Uses the nsold function by C.T. Kelley.
%	Available here: https://github.com/ctkelley/ctkelley.github.io/blob/b0aa502875f9409a16254d81de37a8d9a7ce0cf4/newton/SOLVERS/nsold.m


% For methods which require one nonlinear system
%iterative_scheme = 'newton';
%iterative_scheme = 'simplified newton 1';
%iterative_scheme = 'simplified newton 2'; 
%iterative_scheme = 'simplified newton 3'; 
iterative_scheme = 'kelley';

% For methods which require two nonlinear systems
%iterative_scheme1 = 'newton';
%iterative_scheme1 = 'simplified newton 0'; % Lobatto IIIA-IIIB only
%iterative_scheme1 = 'simplified newton';
%iterative_scheme1 = 'simplified newton 2';
%iterative_scheme1 = 'simplified newton 3';
iterative_scheme1 = 'kelley';
%iterative_scheme2 = 'newton';
%iterative_scheme2 = 'simplified newton';
iterative_scheme2 = 'kelley';

% -------
% Problem
% -------

% pendulum: The plane pendulum problem.  k is linear in u.
% alternate pendulum: A pendulum-spring model with damping.  k is linear in u.
% linear: An index 3 DAE in which k is linear in u.
% nonlinear: An example where k is nonlinear in u.  Requires that
% 	plot_str be set to component.

%problem_str = 'pendulum';
%problem_str = 'alternate pendulum';
problem_str = 'linear';
%problem_str = 'nonlinear';

% -------
% Method
% -------

% Note that there is a separate method_str for Lobatto IIIC and
% projected Lobatto IIIC.

%method_str = 'Gauss';
%method_str = 'Radau IIA';
%method_str = 'projected Radau IIA';
%method_str = 'Lobatto IIIC';
method_str = 'Lobatto IIIA-IIIB';
%method_str = 'projected Lobatto IIIC';
%method_str = 'projected Lobatto III';

% ----------------
% Number of stages
% ----------------

%s=2; % Radau IIA, Lobatto IIIA-IIIB
s=3; % Radau IIA, projected Radau IIA, Gauss, Lobatto IIIC, projected Lobatto IIIC, Lobatto IIIA-IIIB, projected Lobatto III

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
% pendulum: Shows how well the pendulum stays on the manifold 0=g(y) and
% 	hidden manifold 0=(g_y f)(y,z).  Also plots position coordates and
% 	energy as a function of time.  Can only be used if
% 	problem_str='pendulum'.
% component: Plot each of y_j, z_j, and u_j vs. t separately.
% linear: Plot constraint violation for the linear DAE
% empty string: No plots.

%plot_str = 'pendulum'; 
%plot_str = 'component'; 
plot_str = 'linear';
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
N = 115;
% ================================================================

% Get only the number of nonlinear systems from get_values_method
method_result = get_values_method(method_str, s);
number_nonlinear_systems = method_result.number_nonlinear_systems;

% Get only the true solution, if applicable, from get_values_problem
problem_result = get_values_problem(problem_str);
true_solution = problem_result.true_solution;

% Make a structured array of the iterative schemes.  
% One element per nonlinear system to solve.
if number_nonlinear_systems==1
	iterative_scheme_struct.iterative_scheme = iterative_scheme;
elseif number_nonlinear_systems==2
	iterative_scheme_struct.iterative_scheme1 = iterative_scheme1; 
	iterative_scheme_struct.iterative_scheme2 = iterative_scheme2; 
else
	error('Error (init): Invalid value of number_nonlinear_systems')
end

% -------------
% Miscellaneous
% -------------

% Check that either plots or a movie is produced, not both
if movie && strlength(plot_str) > 0
	error('Error (init): movie=true requires that plot_str is empty')
end

% check that the plot and movie options make sense for the problem
if strcmpi(plot_str, 'pendulum') && ~strcmpi(problem_str, 'pendulum')
	error('Error (init): The plotting option pendulum can only be chosen for the pendulum problem.')
elseif strcmpi(problem_str, 'nonlinear') && ~strcmpi(plot_str, 'component')
	error('Error (init): The nonlinear problem requires that plot_str be set to component.')
end

if movie && ~strcmpi(problem_str, 'pendulum')
	error('Error (init): movie=true requires that the problem to be solved is the pendulum')
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
plots_params.iterative_scheme_struct = iterative_scheme_struct;
plots_params.number_nonlinear_systems = number_nonlinear_systems;
plots_params.plot_iterations = plot_iterations;

% append true solution to plots_params, if applicable
if strcmpi(problem_str, 'nonlinear') || strcmpi(problem_str, 'linear')
	plots_params.true_solution = true_solution;
end


