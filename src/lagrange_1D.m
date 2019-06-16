function coeff = lagrange_1D( xc_frng, xc_don, ind, p)
% This function computes the coefficient corresponding to ind for 
% interpolating the solution at the fringe node from the donor nodes in 1D.

% determine coefficient corresponding to ind node
den = 1.0;
num = 1.0;
for jnd=1:p
    if(jnd~=ind)
         den = den * (xc_don(ind) - xc_don(jnd));
         num = num * (    xc_frng - xc_don(jnd));
    end
end
coeff = num/den;

end
