function Khat_YZU = computeKhat(Yhat,Z,Uhat,k,ny,nz,nu,s,y0)
	% function Khat_YZU = computeKhat(Yhat,Z,Uhat,k,ny,nz,nu,s,y0)
    % Returns Khat(Yhat,Z,Uhat) = [k(Y0,Z0,U0); ...;k(Y_{s-1},Z_{s-1},U_{s-1})].
	% For use with the projected Lobatto III method.
    %
%% BOOKMARK
    % Inputs:
    %   Yhat: Column vector [Y_2; ...; Y_{s+1}]
    %   Z: Column vector [Z_1; ...; Z_s]
    %   Uhat: Column vector [U_0; ...; U_{s-1}]
    %   k: RHS of z'=k(y,z,u)
    %   ny: dimension of Y_i 
    %   nz: dimension of Z_i
    %   nu: dimension of U_i
    %   s: Number of stages
    %   y0: Previous step of y-component.  
    %
    % Output:
    %   Khat_YZU: The vector [k(Y_0,Z_0,U_0); ...; k(Y_{s-1},Z_{s-1},U_{s-1})]

	% check that y0 is a column vector
	if size(y0,2) > 1
		error('Error (computeKhat): y0 must be a column vector')
	end

    % Reshape into matrices
    Yhat_mat = reshape(Yhat, ny, s); % the matrix [Y_2 ... Y_{s+1}]
    Z_mat = reshape(Z, nz, s); % the matrix [Z1 ... Zs]
    Uhat_mat = reshape(Uhat, nu, s); % the matrix [U_0 ... U_{s-1}]

	Yhat_mat = Yhat_mat(:,1:s-2); % [Y_2 ... Y_{s-1}]
	Yhat_mat = [y0 y0 Yhat_mat]; % [Y_0 Y_1 ... Y_{s-1}]

	Z0 = Z_mat(:,1); % Z_0 = Z_1
	Z_mat = Z_mat(:,1:s-1); % [Z_1 ... Z_{s-1}]
	Z_mat = [Z0 Z_mat]; % [Z_0 Z_1 ... Z_{s-1}]

    % Convert to cell arrays, one column vector per cell
    Y_cells = mat2cell(Yhat_mat, ny, ones(1, s));
    Z_cells = mat2cell(Z_mat, nz, ones(1, s));
    U_cells = mat2cell(Uhat_mat, nu, ones(1, s));

    % Apply k to each 3-tuple of cells
    Kcells = cellfun(k, Y_cells, Z_cells, U_cells, 'UniformOutput', false);

    % Concatenate into a single column vector
    Khat_YZU = vertcat(Kcells{:});
end


