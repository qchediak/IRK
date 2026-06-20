function G_Y = computeG(Y,g,ny,s)
	% function G_Y = computeG(Y,g,ny,s)
	% returns G(Y), where Y=[Y_1;...;Y_s] and
	% G(Y) = [g(Y_1);...;g(Y_s)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	g: RHS of the constraint equation 0=g(y)
	%	ny: common dimension of Y1,...,Ys
	% 	s: number of internal stages
	%
	% Outputs:
	%	G_Y: This is the column vector 
	% 		G(Y) = [g(Y_1);...;g(Y_s)].

    Y_matrix = reshape(Y, ny, s); % convert Y to a matrix [Y1 Y2 ... Ys]
    Y_cells = mat2cell(Y_matrix, ny, ones(1, s)); % convert matrix to cell array structure
    G_cells = cellfun(g, Y_cells, 'UniformOutput', false); % apply g to each cell
    G_Y = vertcat(G_cells{:}); % concatenate into a single column
end


