function dFdY_Y = compute_dFdY_ODE(Y,dfdy,ny,s)
	%function dFdY_Y = compute_dFdY_ODE(Y,dfdy,ny,s)
	% Returns the Jacobian matrix (dF/dY)(Y), 
	% where Y=[Y1;...;Ys] and F(Y)=[f(Y1);...;f(Ys)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	dfdy: The Jacobian matrix df/dy
	%	ny: common dimension of Y1,...,Ys
	% 	s: number of internal stages
	%
	% Outputs:
	%	dFdY_Y: The (s*ny) x (s*ny) Jacobian matrix (dF/dY)(Y)

	% split Y into  cell arrays of ny x 1 vectors
	Y_cells = mat2cell(Y, ny * ones(s,1), 1);

	% apply df/dy to each cell
	dfdy_cells = cellfun(dfdy, Y_cells, 'UniformOutput', false);

	% Make the block diagonal matrix
	dFdY_Y = blkdiag(dfdy_cells{:});
end



