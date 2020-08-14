function coeff = coeff_rbf(frng_coords,supp_coords,ov_info)
% assemble interpolation coefficients corresponding to gaussian RBF

% determine coefficients based on polynomial support desired
if ov_info('poly order') < 0
    coeff = coeff_rbf_classic(frng_coords,supp_coords,ov_info);
else
    coeff = coeff_rbf_augmented(frng_coords,supp_coords,ov_info);
end

end

%% assemble interpolation coefficients based on calssic RBF - based on Mesh Free methods by Liu

function coeff = coeff_rbf_classic(frng_coords,supp_coords,ov_info)

num_intrp_coeff = size(supp_coords,1); % number of interpolation coefficients

% allocate space for moment matrix
R_Q = zeros(num_intrp_coeff,num_intrp_coeff);

% assemble radial basis for fringe coordinates
R = gaussian_basis(frng_coords,supp_coords,ov_info('shape param'));

for n = 1:num_intrp_coeff
    
    % assemble radial basis momemt matrix
    switch ov_info('intrp shape')
        case "gaussian"
            R_Q(n,:) = (gaussian_basis(supp_coords(n,:),supp_coords,ov_info('shape param')))';
            
        otherwise
            error('Interpolation shape type not supported;check overset info in driver');
    end
end

coeff = (R'/R_Q)'; % evaluate radial basis coefficient array

end

%% assemble interpolation coefficients based on augmented RBF - based on Mesh Free methods by Liu

function coeff = coeff_rbf_augmented(frng_coords,supp_coords,ov_info)

num_intrp_coeff = size(supp_coords,1); % number of interpolation coefficients

% determine coefficients corresponding to desired polynomial support
poly_order = ov_info('poly order'); 
num_poly_coeff = factorial(poly_order+2)/factorial(2)/factorial(poly_order);

% allocate space for moment matrices
R_Q = zeros(num_intrp_coeff,num_intrp_coeff);
P_m = zeros(num_intrp_coeff,num_poly_coeff);

% assemble radial basis for fringe coordinates
R = gaussian_basis(frng_coords,supp_coords,ov_info('shape param'));

% assemble polynomial basis for frng_coords
p = assemble_poly_coeff(frng_coords,num_poly_coeff);

for n = 1:num_intrp_coeff
    
    % assemble radial basis momemt matrix
    switch ov_info('intrp shape')
        case "gaussian"
            R_Q(n,:) = (gaussian_basis(supp_coords(n,:),supp_coords,ov_info('shape param')))';
            
        otherwise
            error('Interpolation shape type not supported;check overset info in driver');
    end

    % assemble polynomial basis moment matrix
    P_m(n,:) = (assemble_poly_coeff(supp_coords(n,:),num_poly_coeff))';
end

P_mT_R_Q_inv = P_m'/R_Q; % pre compute array to be used in S_a and S_b
S_b = (P_mT_R_Q_inv*P_m)\P_mT_R_Q_inv;
S_a = eye(num_intrp_coeff)/R_Q - P_mT_R_Q_inv'*S_b;

coeff = (R'*S_a + p'*S_b)'; % evaluate radial basis coefficient array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Jay's implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M = zeros(num_intrp_coeff+num_poly_coeff);
% B = zeros(num_intrp_coeff+num_poly_coeff,1);
% 
% for i = 1:num_intrp_coeff
%     
%     for j = 1:num_intrp_coeff  
%         dd = dot(supp_coords(i,:)-supp_coords(j,:),supp_coords(i,:)-supp_coords(j,:));
%         M(i,j) = exp(-dd);    
%     end
%     M(i,num_intrp_coeff+1:end) = assemble_poly_coeff(supp_coords(i,:),num_poly_coeff);
%     M(num_intrp_coeff+1:end,i) = M(i,num_intrp_coeff+1:end);
%     
%     dd = dot(frng_coords-supp_coords(i,:),frng_coords-supp_coords(i,:));
%     B(i) = exp(-dd); 
% end
% 
% B(num_intrp_coeff+1:end) = assemble_poly_coeff(frng_coords,num_poly_coeff);
% 
% coeff_test = M\B;
% coeff = coeff_test(1:num_intrp_coeff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% function to assemble vector of polynomial support
function p = assemble_poly_coeff(coord,num_poly_coeff)

x = coord(1); % x coordinate
y = coord(2); % y coordinate

p = zeros(num_poly_coeff,1); % allocate memory

% determine donor grid
switch num_poly_coeff
    case 1 % constant 
        p(1)  = 1;
        
    case 3 % 1st order polynomial
        p(1)  = 1;
        p(2)  = x;
        p(3)  = y;

    case 6 % 2nd order polynomial
        p(1)  = 1;
        p(2)  = x;
        p(3)  = y;
        p(4)  = x*x;
        p(5)  = x*y;
        p(6)  = y*y;
        
    case 10 % 3rd order polynomial
        p(1)  = 1;
        p(2)  = x;
        p(3)  = y;
        p(4)  = x*x;
        p(5)  = x*y;
        p(6)  = y*y;
        p(7)  = x*x*x;
        p(8)  = x*x*y;
        p(9)  = x*y*y;
        p(10) = y*y*y;

    otherwise
        error('Order of polynomial support not yet included;check overset info in driver');
end

end