%% Make cylinders 

% coord = centre of cylinder based on xyz coordinate
% R = radius of cylinder
% L = length of cylinder
% nCS = number of cross sections >> more detailed rings lengthwise 
% nNodes = number of nodes >> more detailed circlular ring 

function cyl_ob = mk_cyl(coord, R, L, nCS, nNodes)

    for c = 1:size(coord,1)
    
        centre = coord(c,:); 

        % Generate cylinder coordinates
        r = R * ones(1, nNodes);
        th = linspace(0, 2*pi, nNodes);
        [x, y] = pol2cart(th, r);

        % Offset coordinates by the center values
        x = x + centre(1);
        y = y + centre(2);
        z = linspace(0, L, nCS)' + (centre(3) - L/2);

        X = repmat(x, nCS, 1);
        Y = repmat(y, nCS, 1);
        Z = repmat(z, 1, nNodes);

        % Generate lid coordinates
        x_lid = zeros(2, nNodes) + centre(1);
        y_lid = zeros(2, nNodes) + centre(2);
        z_lid = repmat([0; L] + (centre(3) - L/2), 1, nNodes);

        % Concatenate lid and cylinder coordinates
        X = [x_lid(1, :); X; x_lid(2, :)];
        Y = [y_lid(1, :); Y; y_lid(2, :)];
        Z = [z_lid(1, :); Z; z_lid(2, :)];

        % Plot the cylinder with lids
        
        dims(1,:,:) = X; 
        dims(2,:,:) = Y; 
        dims(3,:,:) = Z; 


        cyl_ob{c} = dims; 

    end 

end 
