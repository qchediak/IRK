function dFdZ_YZ = compute_dFdZ(Y,Z,dfdz,ny,nz,s)
	% function dFdZ_YZ = compute_dFdZ(Y,Z,dfdz,ny,nz,s)
	% Returns the Jacobian matrix (dF/dZ)(Y,Z), 
	% where Y=[Y1;...;Ys], Z=[Z1; ...; Zs],
	% and F(Y,Z)=[f(Y1,Z1);...;f(Ys,Zs)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	Z: Column vector of internal stages [Z1;...;Zs]
	%	dfdy: The partial derivative df/dy
	%	ny: common dimension of Y1,...,Ys
	%	nz: common dimension of Z1,...,Zs
	% 	s: number of internal stages
	%
	% Outputs:
	%	dFdY_Y: The (s*ny) x (s*nz) matrix (dF/dZ)(Y,Z)

	% split Y and Z into  cell arrays of ny x 1 or nz x 1 vectors
	Y_cells = mat2cell(Y, ny * ones(s,1), 1);
	Z_cells = mat2cell(Z, nz * ones(s,1), 1);

	% apply df/dy to each cell
	dfdz_cells = cellfun(dfdz, Y_cells, Z_cells, 'UniformOutput', false);

	% Make the block diagonal matrix
	dFdZ_YZ = blkdiag(dfdz_cells{:});
end




