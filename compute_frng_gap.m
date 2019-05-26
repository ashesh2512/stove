function compute_frng_gap(mesh1,mesh2,nb_ib1,bg_ib1,nb_ib2,bg_ib2)
% compute distance between closest fringes in near body and background
% meshes

coords1 = mesh1{2,1}; % coordinates for mesh 1
coords2 = mesh2{2,1}; % coordinates for mesh 2

% skip if mesh 2 not treated as near body
if ~isempty(coords2((nb_ib2==-1)))
    % distance between innermost fringes along bottom side with mesh 2 as near
    % body
    inner_frng_top_plane2 = nb_ib2==-1 ...
                          & coords2(:,2)< min(coords1(bg_ib2==-1,2)) ...
                          & coords2(:,1)>=min(coords1(bg_ib2==-1,1)) ...
                          & coords2(:,1)<=max(coords1(bg_ib2==-1,1));
    frng_gap_bot2 = abs(max(coords2(inner_frng_top_plane2,2)) - min(coords1(bg_ib2==-1,2)));
    
    % distance between innermost fringes along right side with mesh 2 as near
    % body
    inner_frng_right_plane2 = nb_ib2==-1 ...
                            & coords2(:,1)> max(coords1(bg_ib2==-1,1)) ...
                            & coords2(:,2)>=min(coords1(bg_ib2==-1,2)) ...
                            & coords2(:,2)<=max(coords1(bg_ib2==-1,2));
    frng_gap_right2 = abs(min(coords2(inner_frng_right_plane2,1)) - max(coords1(bg_ib2==-1,1)));
    
    % distance between innermost fringes along top side with mesh 2 as near
    % body
    inner_frng_top_plane2 = nb_ib2==-1 ...
                          & coords2(:,2)> max(coords1(bg_ib2==-1,2)) ...
                          & coords2(:,1)>=min(coords1(bg_ib2==-1,1)) ...
                          & coords2(:,1)<=max(coords1(bg_ib2==-1,1));
    frng_gap_top2 = abs(min(coords2(inner_frng_top_plane2,2)) - max(coords1(bg_ib2==-1,2)));
    
    % distance between innermost fringes along left side with mesh 2 as near
    % body
    inner_frng_left_plane2 = nb_ib2==-1 ...
                           & coords2(:,1)< min(coords1(bg_ib2==-1,1)) ...
                           & coords2(:,2)>=min(coords1(bg_ib2==-1,2)) ...
                           & coords2(:,2)<=max(coords1(bg_ib2==-1,2));
    frng_gap_left2 = abs(max(coords2(inner_frng_left_plane2,1)) - min(coords1(bg_ib2==-1,1)));
    
    fprintf('\n near body mesh 2 fringe gap: bottom: %e, right: %e, top: %e, left:%e', ...
            frng_gap_bot2, frng_gap_right2, frng_gap_top2, frng_gap_left2);
    fprintf('\n');
end

% return if mesh 1 not treated as near body
if isempty(coords1((nb_ib1==-1)))
    return
end

% distance between innermost fringes along bottom side with mesh 1 as near
% body
inner_frng_bot_plane1 = nb_ib1==-1 ...
                      & coords1(:,2)< min(coords2(bg_ib1==-1,2)) ...
                      & coords1(:,1)>=min(coords2(bg_ib1==-1,1)) ...
                      & coords1(:,1)<=max(coords2(bg_ib1==-1,1));       
frng_gap_bot1 = abs(max(coords1(inner_frng_bot_plane1,2)) - min(coords2(bg_ib1==-1,2)));

% distance between innermost fringes along right side with mesh 1 as near
% body
inner_frng_right_plane1 = nb_ib1==-1 ...
                        & coords1(:,1)> max(coords2(bg_ib1==-1,1)) ...
                        & coords1(:,2)>=min(coords2(bg_ib1==-1,2)) ...
                        & coords1(:,2)<=max(coords2(bg_ib1==-1,2));                        
frng_gap_right1 = abs(min(coords1(inner_frng_right_plane1,1)) - max(coords2(bg_ib1==-1,1)));

% distance between innermost fringes along top side with mesh 1 as near
% body
inner_frng_top_plane1 = nb_ib1==-1 ...
                      & coords1(:,2)> max(coords2(bg_ib1==-1,2)) ...
                      & coords1(:,1)>=min(coords2(bg_ib1==-1,1)) ...
                      & coords1(:,1)<=max(coords2(bg_ib1==-1,1));       
frng_gap_top1 = abs(min(coords1(inner_frng_top_plane1,2)) - max(coords2(bg_ib1==-1,2)));

% distance between innermost fringes along left side with mesh 1 as near
% body
inner_frng_left_plane1 = nb_ib1==-1 ...
                       & coords1(:,1)< min(coords2(bg_ib1==-1,1)) ...
                       & coords1(:,2)>=min(coords2(bg_ib1==-1,2)) ...
                       & coords1(:,2)<=max(coords2(bg_ib1==-1,2));                        
frng_gap_left1 = abs(max(coords1(inner_frng_left_plane1,1)) - min(coords2(bg_ib1==-1,1)));

% print the fringe gaps
fprintf('\n near body mesh 1 fringe gap: bottom: %e, right: %e, top: %e, left:%e', ...
        frng_gap_bot1, frng_gap_right1, frng_gap_top1, frng_gap_left1);
fprintf('\n');

end