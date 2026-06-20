% 2026 Jan 24
% Find orders of starting approximations for ODEs

clear; close all

path(pathdef); % Reset search path to saved default
addpath('../shared_functions'); % add the shared functions folder

init

% -----------------------------------------------------
% Initialize the error vectors as vectors of NaN values
% -----------------------------------------------------

% If error_plotting_Y='composite', then Y0_errs(i)
% is the error in Y_{n+1}^0 at the stepsize h=stepsizes(i).
%
% So, Y0_errs is the column vector of errors in Y_{n+1}^0
% at each stepsize.
%
% If error_plotting_Y='component', then Y0_errs(i,j) is the error of
% Y_{n+1,j}^0 at stepsize h=stepsizes(i). 
%
% So, Y0_errs(:,j) is the column vector of errors of Y_{n+1,j}^0
% at all stepsizes, and Y0_errs(i,:) is the row vector of all
% the errors at the fixed stepsize h=stepsizes(i).

if strcmpi(error_plotting_Y, 'composite')
	Y0_errs = NaN(stepsizes_len,1);
elseif strcmpi(error_plotting_Y, 'component') 
	% Note for Lobatto IIIA-IIIB Y0_errs(:,1) will remain NaN
	Y0_errs = NaN(stepsizes_len,s);
else
	error('Error (main): Invalid value of error_plotting_Y')
end

% ----------
% Main loop
% ----------

for i=1:stepsizes_len

	h = stepsizes(i);

	% The main function call
	result = find_error_SA_ODE(h, params, find_error_SA_option);

	% unpack result
	steps = result.steps;
	err_vecs = result.err_vecs;

	% unpack steps
	yns = steps.yns;

	% unpack err_vecs
	Y0_err_vec = err_vecs.Y0_err_vec;

	% ---------------------
	% update error vectors
	% ---------------------

	% Find errors for Y0
	%
	% Keep errors in either component or composite form
	% depending on error plotting option set in init file.
	% If the option is set to 'composite', take the norm of the entire
	% error vector.  Otherwise, take the norm of each individual
	% component.

	% Y component
	if strcmpi(error_plotting_Y, 'composite')
		Y0_errs(i,1) = norm(Y0_err_vec);
	elseif strcmpi(error_plotting_Y, 'component') 
		for j=1:s % j is stage number, i corresponds to value of h
			Y0_errs(i,j) = norm(Y0_err_vec((j-1)*ny+1:j*ny));
		end
	else 
		error('Error (main): Invalid value of error_plotting_Y')
	end
end

% make errs_struct
errs_struct.Y0_errs = Y0_errs;

% print info
print_starting_approx_info_ODE(stepsizes, print_plot_params)

% plot the errors
result_errs = plot_all_errors_ODE(stepsizes, errs_struct, print_plot_params);


