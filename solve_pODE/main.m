% 2025 Dec 14
% main script for solve_pODE

clear; close all

path(pathdef); % Reset search path to saved default
addpath('../shared_functions'); % add the shared functions folder

init 

% the main function call
result = find_trajectory_pODE(iterative_scheme, problem_str, method_str, starting_approx_str, h, s, tol, N);

% unpack result
yns = result.yns;
zns = result.zns;
iter_vec = result.iter_vec;

if strcmpi(problem_str, 'pendulum')
	E = result.E; % energy
end

% Print info.
% Note that print_info_ODE can be used for both ODEs and partitioned ODEs.
print_info_ODE(method_str, problem_str, s, h, N, iter_vec, iterative_scheme, starting_approx_str, tol)

% ---------------
% Plots and Movie
% ---------------

% append energy values to plots_params, if applicable
if strcmpi(problem_str, 'pendulum')
	plots_params.E = E;
end

% pendulum only
% Note that make_movie can be used for both partitioned ODEs and index 3
% DAEs.
if movie 
	make_movie(yns,zns,plots_params)
end

% Plots
if strlength(plot_str) > 0 
	make_plots_pODE(yns,zns,iter_vec,plots_params)
end

