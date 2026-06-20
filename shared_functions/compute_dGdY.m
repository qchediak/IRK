function dGdY_Y = compute_dGdY(Y,dgdy,ny,s)
	% function dGdY_Y = compute_dGdY(Y,dgdy,ny,s)
	% Returns the Jacobian matrix (dG/dY)(Y), where Y=[Y_1;...;Y_s], and
	% G(Y)=[g(Y1);...;g(Ys)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	ny: common dimension of Y1,...,Ys
	%	dgdy: the derivative dg/dy 
	%		(i.e., the derivative of g with respect to its argument)
	% 	s: number of internal stages
	%
	% Outputs:
	%	dGdY_Y: The (nu*s)x(ny*s) matrix (dG/dY)(Y)

	% split Y into a cell array of ny x 1 vectors
	Y_cells = mat2cell(Y, ny * ones(s,1), 1);

	% apply transpose and scaling to each block
	dgdy_cells = cellfun(dgdy, Y_cells, 'UniformOutput', false);

	% Make the block diagonal matrix
	dGdY_Y = blkdiag(dgdy_cells{:});
end


