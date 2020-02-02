function [fdof,fld_dof,dnd,pnd,nd_dof_map_pnd] = get_dof_status(mesh,bc,nd_dof_map,ov_info,iblank)
% apply initial and boundary conditions based on problem selected
%
% Input:  mesh       - mesh object created uing meshgen
%         bc         - boundary conditions
%         nd_dof_map - map between node and dofs
%
% Output: fdof     - array of non-constrained degrees of freedom
%         fld_dof  - degrees of freedom corresponding to field nodes
%         dnd      - array of dirichlet nodes
%         pnd      - array of periodic nodes (only slave constrained)

% extract edge topologies
edge_topo_x = mesh{2,2};
edge_topo_y = mesh{2,4};

% extract sidesets
bottom_ss = mesh{2,7};
right_ss  = mesh{2,8};
top_ss    = mesh{2,9};
left_ss   = mesh{2,10};

% collect all cbc nodes
dnd = []; pnd_mas = []; pnd_slv = [];

if (bc('bottom') == "dirichlet")
    dnd = [dnd; reshape(edge_topo_x(bottom_ss,:),[],1)];
elseif (bc('bottom') == "per_slave")
    pnd_slv = [pnd_slv; unique(reshape(edge_topo_x(bottom_ss,:),[],1))];
    pnd_mas = [pnd_mas; unique(reshape(edge_topo_x(top_ss,   :),[],1))];
end
if (bc('right') == "dirichlet")
    dnd = [dnd; reshape(edge_topo_y(right_ss,:),[],1)];
elseif (bc('right') == "per_slave")
    pnd_slv = [pnd_slv; unique(reshape(edge_topo_y(right_ss,:),[],1))];
    pnd_mas = [pnd_mas; unique(reshape(edge_topo_y(left_ss, :),[],1))];
end
if (bc('top') == "dirichlet")
    dnd = [dnd; reshape(edge_topo_x(top_ss,:),[],1)];
elseif (bc('top') == "per_slave")
    pnd_slv = [pnd_slv; unique(reshape(edge_topo_x(top_ss,    :),[],1))];
    pnd_mas = [pnd_mas; unique(reshape(edge_topo_x(bottom_ss, :),[],1))];
end
if (bc('left') == "dirichlet")
    dnd = [dnd; reshape(edge_topo_y(left_ss,:),[],1)];
elseif (bc('left') == "per_slave")
    pnd_slv = [pnd_slv; unique(reshape(edge_topo_y(left_ss ,:),[],1))];  
    pnd_mas = [pnd_mas; unique(reshape(edge_topo_y(right_ss,:),[],1))];   
end

% if all edges are periodic, the diagonal nodes
% (nodes common to the master boundaries and slave boundaries) are paired too! 
if (bc('bottom') == "per_slave" && bc('left') == "per_slave")
    pnd_mas = [pnd_mas; size(nd_dof_map,1)];
    pnd_slv = [pnd_slv; 1];
elseif (bc('bottom') == "per_slave" && bc('right') == "per_slave")
    mas_nd = unique(reshape(edge_topo_y(left_ss  ,:),[],1));
    pnd_mas = [pnd_mas; mas_nd(end)];
    slv_nd = unique(reshape(edge_topo_x(bottom_ss,:),[],1)); slv_nd = slv_nd(end);
    pnd_slv = [pnd_slv; slv_nd(end)];
elseif (bc('top') == "per_slave" && bc('left') == "per_slave")
    slv_nd = unique(reshape(edge_topo_y(left_ss  ,:),[],1));
    pnd_slv = [pnd_slv; slv_nd(end)];
    mas_nd = unique(reshape(edge_topo_x(bottom_ss,:),[],1));
    pnd_mas = [pnd_mas; mas_nd(end)];
elseif (bc('top') == "per_slave" && bc('right') == "per_slave")
    pnd_mas = [pnd_mas; 1];
    pnd_slv = [pnd_slv; size(nd_dof_map,1)];
end

dnd = unique(dnd);
pnd = [pnd_mas pnd_slv];

% ensure a slave node isn't also a master node
% multiple slaves can contribute to the same master but one slave cannot
% contribute to multiple masters
nd_dof_map_pnd = nd_dof_map; % initialize modified node to dof map incorporating periodic bc

if(~isempty(pnd))
    [~,ind] = intersect(pnd(:,1), pnd(:,2));
    pnd(ind,:) = [];
    if(sort(pnd(:,2)) ~= unique(pnd(:,2)))
        error('Slave donating to multiple masters. Fix periodic pairing.');
    end
    nd_dof_map_pnd(pnd(:,2),:) = nd_dof_map(pnd(:,1),:);
end

% extract array of free dofs
fdof = unique(setdiff(nd_dof_map, nd_dof_map([dnd;pnd_slv])));

% set array of field dofs
fld_dof = nd_dof_map(iblank==1,:);

% for a fully decoupled overset solve constrain overset nodes too
if ~(ov_info('solve type') == "decoupled" && ov_info('fringe update') == "direct")
    return
end

% update free dof array to exclude dofs on fringe nodes
fdof = setdiff(fdof, nd_dof_map(iblank==-1,:));

end
