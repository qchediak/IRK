function dKdU_YZU = compute_dKdU(Y,Z,U,dkdu,ny,nz,nu,s)
	% function dKdU_YZU = compute_dKdU(Y,Z,U,dkdu,ny,nz,nu,s)
	%
	% Returns the Jacobian matrix (dK/dU)(Y,Z,U), 
	% where Y=[Y_1;...;Y_s], Z=[Z_1;...;Z_s], U=[U1;...;Us], 
	% and K(Y,Z,U)=[k(Y1,Z1,U1);...;k(Ys,Zs,Us)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	Z: Column vector of internal stages [Z1;...;Zs]
	%	U: Column vector of internal stages [U1;...;Us]
	%	dkdu: partial deriviative dk/du
	%	ny: common dimension of Y1,...,Ys, 
	%	nz: common dimension of Z1,...,Zs, 
	%	nu: common dimension of U1,...,Us
	% 	s: number of internal stages
	%
	% Outputs:
	%	dKdU_YZU: The (s*nz)x(s*nu) matrix (dK/dU)(Y,Z,U).

	% split Y, Z, U into  cell arrays of nyx1, nzx1, or nux1 vectors
	Y_cells = mat2cell(Y, ny * ones(s,1), 1);
	Z_cells = mat2cell(Z, nz * ones(s,1), 1);
	U_cells = mat2cell(U, nu * ones(s,1), 1);

	% apply dk/du to each cell
	dkdu_cells = cellfun(dkdu, Y_cells, Z_cells, U_cells, 'UniformOutput', false);

	% Make the block diagonal matrix
	dKdU_YZU = blkdiag(dkdu_cells{:});
end


