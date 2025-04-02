%% Scatter 3D Lazy 

% use this by putting a 3d matrix in to save time writing 


function scatter3_lazy(matrix,c) 


if ~exist('c','var')
    % third parameter does not exist, so default it to something
    c = 'k';
end

scatter3(matrix(:,1), matrix(:,2), matrix(:,3), c)

end