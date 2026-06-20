function Fhat_YZ = computeFhat(Yhat,Z,f,ny,nz,s,y0)
	% function F_YZ = computeF(Yhat,Z,f,ny,nz,s,y0)
	% returns Fhat(Yhat,Z), where Yhat=[Y_2;...;Y_{s+1}],
	% Z=[Z_1;...;Z_s], and Fhat(Yhat,Z) = [f(Y_1,Z_1);...;f(Y_s,Z_s)].
	% For use with the projected Lobatto III method.
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y_2;...;Y_{s+1}]
	%	Z: Column vector of internal stages [Z_1;...;Z_s]
	%	f: RHS of y'=f(y,z)
	%	ny: common dimension of Y_2,...,Y_{s+1}
	%	nz: common dimension of Z_1,...,Z_s
	% 	s: number of internal stages
	% 	y0: Equal to Y1; the previous step.  Considered a parameter
	% 		here.
	%
	% Outputs:
	%	Fhat_YZ: This is the column vector 
	% 		Fhat(Yhat,Z) = [f(Y_1,Z_1);...;f(Y_s,Z_s)].

	% check that y0 is a column vector
	if size(y0,2) > 1
		error('Error (computeFhat): y0 must be a column vector')
	end

    Yhat_matrix = reshape(Yhat, ny, s); % convert Y to a matrix [Y_2 Y_3 ...  Y_{s+1}]
    Z_matrix = reshape(Z, nz, s); % convert Z to a matrix [Z1 Z2 ... Zs]

	Yhat_matrix = Yhat_matrix(:,1:s-1); % only keep [Y2 ... Ys]
	Yhat_matrix = [y0 Yhat_matrix]; % now it's [y0 Y2 ... Ys]

    Yhat_cells = mat2cell(Yhat_matrix, ny, ones(1, s)); % convert matrix to cell array structure
    Z_cells = mat2cell(Z_matrix, nz, ones(1, s)); % convert matrix to cell array structure
    Fhat_cells = cellfun(f, Yhat_cells, Z_cells, 'UniformOutput', false); % apply f to each pair of cells
    Fhat_YZ = vertcat(Fhat_cells{:}); % concatenate into a single column
end


