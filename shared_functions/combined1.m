function [f,jac] = combined1(x, H, DH)
	% function [f,jac] = combined1(x, H, DH)
	% for use with nsold
    f = H(x);

    if nargout == 2
        jac = DH(x);
    end
end
