function coeff = compute_frg_coeff(frng_coords,donor_nd_coords,ov_info)
% assemble interpolation coefficients corresponding to fringe nodes

% determine contribution based on donor grid type
switch ov_info('donor grid')
    case "tensor"
        coeff = coeff_tensor(frng_coords,donor_nd_coords,ov_info('intrp order'));
    
    case "radial"
        coeff = coeff_radial(frng_coords,donor_nd_coords,ov_info);

    otherwise
        error('Do not recognize donor grid type; check overset info in driver');
end

end

%% assemble interpolation coefficients based on tensor grid

function coeff = coeff_tensor(frng_coords,donor_nd_coords,p)

coeff = zeros(size(donor_nd_coords,1),1); % allocate memory for coefficients

xcoords = unique(donor_nd_coords(:,1)); % determine 1D unique x coordinates
ycoords = unique(donor_nd_coords(:,2)); % determine 1D unique y coordinates

intrp_order = p+1; % determine interpolation order

for iy = 1:intrp_order
    for ix = 1:intrp_order
        % determine donor index in 1D
        donor_ind = (iy-1)*intrp_order + ix;

        % compute coefficient based on lagarange interpolation
        coeff(donor_ind,1) = ( lagrange_1D(frng_coords(1),xcoords,ix,intrp_order)  ) ...
                           * ( lagrange_1D(frng_coords(2),ycoords,iy,intrp_order)  );
    end
end

end

%% assemble interpolation coefficients based on radial grid

function coeff = coeff_radial(frng_coords,donor_nd_coords,ov_info)

% determine interpolation type
switch ov_info('intrp type')
    case "rbf"
        coeff = coeff_rbf(frng_coords,donor_nd_coords,ov_info);

    otherwise
        error('Do not recognize dinterpolation type; check overset info in driver');
end

end