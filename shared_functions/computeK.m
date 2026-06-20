function K_YZU = computeK(Y,Z,U,k,ny,nz,nu,s)
	% function K_YZU = computeK(Y,Z,U,k,ny,nz,nu,s)
    % returns K(Y,Z,U) = [k(Y1,Z1,U1); ...;k(Ys,Zs,Us)]
    %
    % Inputs:
    %   Y: Column vector [Y1; Y2; ...; Ys]
    %   Z: Column vector [Z1; Z2; ...; Zs]
    %   U: Column vector [U1; U2; ...; Us]
    %   k: RHS of z'=k(y,z,u)
    %   ny: dimension of Y_i 
    %   nz: dimension of Z_i
    %   nu: dimension of U_i
    %   s: Number of stages
    %
    % Output:
    %   K_YZU: The column vector [k(Y1,Z1,U1); ...; k(Ys,Zs,Us)]

    % Reshape into matrices
    Ymat = reshape(Y, ny, s); % convert Y to a matrix [Y1 Y2 ... Ys]
    Zmat = reshape(Z, nz, s); % convert Z to a matrix [Z1 Z2 ... Zs]
    Umat = reshape(U, nu, s); % convert U to a matrix [U1 U2 ... Us]

    % Convert to cell arrays, one column vector per cell
    Y_cells = mat2cell(Ymat, ny, ones(1, s));
    Z_cells = mat2cell(Zmat, nz, ones(1, s));
    U_cells = mat2cell(Umat, nu, ones(1, s));

    % Apply k to each triple of cells
    Kcells = cellfun(k, Y_cells, Z_cells, U_cells, 'UniformOutput', false);

    % Concatenate into a single column vector
    K_YZU = vertcat(Kcells{:});
end

