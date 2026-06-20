% 2025 Dec 14
% main script for solve_index3

clear; close all

path(pathdef); % Reset search path to saved default
addpath('../shared_functions'); % add the shared functions folder

init 

% the main function call
result = find_trajectory(iterative_scheme_struct, problem_str, method_str, starting_approx_str, h, s, tol, N);

% unpack result
yns = result.yns;
zns = result.zns;
uns = result.uns; 
iterations_struct = result.iterations_struct;

if strcmpi(problem_str, 'pendulum')
	E = result.E; % energy
end

% Print info
print_info(method_str, problem_str, s, h, N, iterations_struct, iterative_scheme_struct, number_nonlinear_systems, starting_approx_str, tol)

% ---------------
% Plots and Movie
% ---------------

% append energy values to plots_params, if applicable
if strcmpi(problem_str, 'pendulum')
	plots_params.E = E;
end

% pendulum only
if movie 
	make_movie(yns,zns,plots_params)
end

% Plots
if strlength(plot_str) > 0 
	make_plots(yns,zns,uns,iterations_struct,plots_params)
end

