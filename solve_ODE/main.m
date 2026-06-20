% 2026 Jan 16
% main file for solve_ODE

clear; close all

path(pathdef); % Reset search path to saved default
addpath('../shared_functions'); % add the shared functions folder

init 

% the main function call
result = find_trajectory_ODE(iterative_scheme, problem_str, method_str, starting_approx_str, h, s, tol, N);

% unpack result
yns = result.yns;
iter_vec = result.iter_vec;

if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1') || strcmpi(problem_str, 'pendulum version 2')
	E = result.E; % energy
end

% Print info
print_info_ODE(method_str, problem_str, s, h, N, iter_vec, iterative_scheme, starting_approx_str, tol)

% ---------------
% Plots and Movie
% ---------------

% append energy values to plots_params, if applicable
if strcmpi(problem_str, 'pendulum') || strcmpi(problem_str, 'pendulum version 1') || strcmpi(problem_str, 'pendulum version 2')
	plots_params.E = E;
end

% pendulum only
if movie 
	make_movie_ODE(yns,plots_params)
end

% Plots
if strlength(plot_str) > 0 
	make_plots_ODE(yns,iter_vec,plots_params)
end

