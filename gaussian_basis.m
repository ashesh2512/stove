function basis = gaussian_basis(center,supp,c)
% function to compute gaussian basis given a center and support
% coordiantes. Sharpness of basis depend on the value of the shape
% parameters, c


dist = vecnorm(center - supp,2,2); % compute eucledian distance

basis = exp(-c.*dist); % gaussian basis

end