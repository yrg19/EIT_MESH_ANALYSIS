
function [cyls_n, ids] = define_elec_nodes(cyls2, Nodes, bound, fig, d)

% This function calculates which nodes in the mesh are within each
% electrode so that they can be assigned the correct value when creating a
% Jacobia. Written 12/02/2024 Y R Gal Shohet 

% cyls2 is a cell array containg the mm coordinates of the cylinders
% ids is a cell array of the indices of the nodes within the electrode. 
% calculated in create_elec_mask 
% Nodes are the nodes we get from mesher (vertices) 
% fig == 1 if you want to plot the locations 

% for i = 1:length(bound)
%     surf_centres(i,:) = mean(Nodes(bound(i,:),:)); % get the mean node coordinate across nodes that that is on a set of boundaries.
% end

surf_centres = Nodes(bound(:),:); 
id_b = bound(:); 
id_taken = []; 


for c= 1:size(cyls2,2)


    cv = reshape(cyls2{c}, 3, []);

    cyl_nodes = [];

    cv = cv / 1000;

    mid = cv(:,end) - cv(:,1); 
    c2 = cv(:,1) + (0.5 * mid); 

    cv = [cv, c2];  
    cylinder_vertices = [cv(:,1), c2, cv(:,end)]; 

    tic
    
    for i = 1:size(cylinder_vertices,2)

        id = []; 

        b = cylinder_vertices(:,i);

        met = surf_centres - b';

        dist_met = vecnorm(met, 2, 2);

        id = find(dist_met < d);

        ign = ismember(id, id_taken); 

        id = id(~ign); 
        
        cyl_nodes = [cyl_nodes; unique(id_b(id))];

        id_taken = [id_taken; id]; 
       
    end
    toc

    cyls_n{c} = Nodes(cyl_nodes,:);

%    scatter3(squeeze(cyls_n{c}(:,1)), squeeze(cyls_n{c}(:,2)), squeeze(cyls_n{c}(:,3)))

    ids{c} = cyl_nodes; 

    disp(['Done calculating cylinder number ', num2str(c)])

    if fig == 1

        scatter3(cyls_n{c}(:,1),cyls_n{c}(:,2), cyls_n{c}(: ,3))
    
        hold on
    end 

end

end