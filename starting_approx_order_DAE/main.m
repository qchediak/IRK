% 2025 Nov 24
% Find orders of starting approximations for index 3 DAEs.

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
% If error_plotting_Y='compoent', then Y0_errs(i,j) is the error of
% Y_{n+1,j}^0 at stepsize h=stepsizes(i). 
%
% So, Y0_errs(:,j) is the column vector of errors of Y_{n+1,j}^0
% at all stepsizes, and Y0_errs(i,:) is the row vector of all
% the errors at the fixed stepsize h=stepsizes(i).
%
% Z0_errs and U0_errs are similarly defined.

if error_plotting_option==1
	if strcmpi(error_plotting_Y, 'composite')
		Y0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_Y, 'component') 
		% Note for Lobatto IIIA-IIIB Y0_errs(:,1) will remain NaN
		Y0_errs = NaN(stepsizes_len,s);
	else
		error('Error (main): Invalid value of error_plotting_Y')
	end

	if strcmpi(error_plotting_Z, 'composite')
		Z0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_Z, 'component')
		Z0_errs = NaN(stepsizes_len,s);
	else
		error('Error (main): Invalid value of error_plotting_Z')
	end

	if strcmpi(error_plotting_U, 'composite')
		U0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_U, 'component') 
		U0_errs = NaN(stepsizes_len,s);
	else
		error('Error (main): Invalid value of error_plotting_U')
	end

	if number_nonlinear_systems==2
		znp0_errs = NaN(stepsizes_len,1);
		if strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
			Unpsp0_errs = NaN(stepsizes_len,1);
		end
	end
elseif error_plotting_option==2 && s==3
	% Note that these are vectors, not cell arrays
	Y10_errs = NaN(stepsizes_len,1);
	Y20_errs = NaN(stepsizes_len,1);
	Y30_errs = NaN(stepsizes_len,1);

	Z0_errs = NaN(stepsizes_len,1);
	U0_errs = NaN(stepsizes_len,1);
else
	error('Error (main): Invalid combination of error_plotting_option and s.')
end

% ----------
% Main loop
% ----------

for i=1:stepsizes_len

	h = stepsizes(i);

	% The main function call
	result = find_error_SA(h, params, find_error_SA_option);

	% unpack result
	steps = result.steps;
	err_vecs = result.err_vecs;

	% unpack steps
	yns = steps.yns;
	zns = steps.zns;
	uns = steps.uns;

	% unpack err_vecs
	Y0_err_vec = err_vecs.Y0_err_vec;
	Z0_err_vec = err_vecs.Z0_err_vec;
	U0_err_vec = err_vecs.U0_err_vec;

	if strcmpi(method_str, 'Lobatto IIIA-IIIB')
		znp0_err_vec = err_vecs.znp0_err_vec;
		Unps0_err_vec = err_vecs.Unps0_err_vec;
	elseif strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
		znp0_err_vec = err_vecs.znp0_err_vec;
		Unpsp0_err_vec = err_vecs.Unpsp0_err_vec;
	end

	% ---------------------
	% update error vectors
	% ---------------------

	if error_plotting_option==1

		% Find errors for Y0, Z0, and U0 (first nonlinear system)
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
			if number_nonlinear_systems==1 || strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
				for j=1:s % j is stage number, i corresponds to value of h
					Y0_errs(i,j) = norm(Y0_err_vec((j-1)*ny+1:j*ny));
				end
			elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
				for j=2:s % j is stage number, i corresponds to value of h
					Y0_errs(i,j) = norm(Y0_err_vec((j-1)*ny+1:j*ny));
				end
			else
				error('Error (main): Invalid method_str')
			end
		else 
			error('Error (main): Invalid value of error_plotting_Y')
		end

		% Z component
		if strcmpi(error_plotting_Z, 'composite')
			Z0_errs(i,1) = norm(Z0_err_vec);
		elseif strcmpi(error_plotting_Z, 'component') 
			for j=1:s % j is stage number, i corresponds to value of h
				Z0_errs(i,j) = norm(Z0_err_vec((j-1)*nz+1:j*nz));
			end
		else 
			error('Error (main): Invalid value of error_plotting_Z')
		end

		% U component
		if strcmpi(error_plotting_U, 'composite')
			if strcmpi(method_str, 'Lobatto IIIA-IIIB')
				U0_err_vec_full = [U0_err_vec; Unps0_err_vec];
				U0_errs(i,1) = norm(U0_err_vec_full);
			else
				U0_errs(i,1) = norm(U0_err_vec);
			end
		elseif strcmpi(error_plotting_U, 'component') 
			if number_nonlinear_systems==1 || strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
				for j=1:s % j is stage number, i corresponds to value of h
					U0_errs(i,j) = norm(U0_err_vec((j-1)*nu+1:j*nu));
				end
			elseif strcmpi(method_str, 'Lobatto IIIA-IIIB')
				for j=1:s-1 % j is stage number, i corresponds to value of h
					U0_errs(i,j) = norm(U0_err_vec((j-1)*nu+1:j*nu));
				end
				U0_errs(i,s) = norm(Unps0_err_vec);
			else
				error('Error (main): Invalid method_str')
			end
		else 
			error('Error (main): Invalid value of error_plotting_U')
		end

		% update znp0 and Unpsp0 error vectors, if applicable
		if number_nonlinear_systems==2
			% Find vector of errors for znp0
			znp0_errs(i) = norm(znp0_err_vec);

			if strcmpi(method_str, 'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
				% Find vector of errors for Unpsp0
				Unpsp0_errs(i) = norm(Unpsp0_err_vec);
			end
		end

	elseif error_plotting_option==2 && strcmpi(method_str, 'projected Lobatto IIIC')

		for j=1:s % j is stage number, i corresponds to value of h
			Y0_errs(i,j) = norm(Y0_err_vec((j-1)*ny+1:j*ny));
		end

		% Z component: Combine Z0 with znp0 
		Z0_errs(i) = norm([Z0_err_vec; znp0_err_vec]);

		% U component: Combine Unp0 with Unpsp0 with U0
		U0_errs(i) = norm([U0_err_vec; Unpsp0_err_vec]);

	else
		error('Error (main): Invalid combination of error_plotting_option and method_str')
	end
end

% make errs_struct
errs_struct.Y0_errs = Y0_errs;
errs_struct.Z0_errs = Z0_errs;
errs_struct.U0_errs = U0_errs;
if number_nonlinear_systems==2 && error_plotting_option==1
	errs_struct.znp0_errs = znp0_errs;
	if strcmpi(method_str,'projected Lobatto IIIC') || strcmpi(method_str, 'projected Radau IIA')
		errs_struct.Unpsp0_errs = Unpsp0_errs;
	end
end

% print info
print_starting_approx_info(stepsizes, print_plot_params)

% plot the errors
result_errs = plot_all_errors(stepsizes, errs_struct, print_plot_params);

