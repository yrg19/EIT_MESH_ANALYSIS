%% Make shank masks

% This script will create a 3D binary mask that represents the shanks in
% the same dimensions as the segments from Seg3D. Created: 19/12/23, Y R
% Gal-Shohet

%% Clear

close all
clear
clc

%% Preamble

% electrode coordinate locations + shank mask output path

inp_file_path = ('/home/yuval/Documents/Creating head mesh/Comsol electrodes/electrode_locations');
out_file_path = ('/home/yuval/Documents/Creating head mesh/e05');

pd = [512, 512, 512];

ps = [0.499022, 0.499022, 0.428571];

res_factor = 5;

pixel_dimensions = pd .* res_factor;

pixel_size = ps ./ res_factor;

addpath('/home/yuval/Packages_and_tools/Random functions')
cd(inp_file_path)

%% Create shanks from data

elec_coord_mm = readmatrix('electrodes_mm_orig.csv');

R = 1.1/ 2; % radius
L = 2.4 ; % length
nCS = 10;
nNodes = 50;

cyls = mk_cyl(elec_coord_mm, R, L, nCS, nNodes);

%% Plot all electrodes

figure

for c = 1 : 64

    surf(squeeze(cyls{c}(1,:,:)), squeeze(cyls{c}(2,:,:)), squeeze(cyls{c}(3,:,:)))
    hold on
end

%% Rotate cylinders to create a shank

% electrode order

rod{1}=(1:6)';
rod{2}=[12, 11, 10, 9, 8, 7, 36, 35, 34, 33]';
rod{3}=[29:31, 54:57]';
rod{4}=[40,41, 17:19, 42]';
rod{5}=[51, 26:28, 52, 53]';
rod{6}=[13:16, 37:39]';
rod{7}=(62:64)';
rod{8}=[24,25,58:61]';
rod{9}=[43, 20, 21, 44:46]';
rod{10}=[22, 47, 23, 48:50]';
rod{11}=32;

% calculate rotations

rot_elecs = computeRotations(rod, L, elec_coord_mm);

% plot rods

figure
for l = 1:length(rod)-1
    subplot(3,4,l)
    shank = elec_coord_mm(rod{l},:)';
    fnplt(cscvn(shank(:,1:end)))

end


%% rotate cylinders

figure

c = 0;
for shankIdx = 1:length(rod) % Iterate over each shank
    for cylIdx = 1:length(rod{shankIdx}) % Iterate over each cylinder in the shank
        electrodeNum = rod{shankIdx}(cylIdx);
        cylVertices = cyls{electrodeNum}; % 3 x 12 x 50 matrix

        % Get rotation angles
        center = elec_coord_mm(electrodeNum,:);
        angles = rot_elecs{shankIdx}(cylIdx, :);
        [rotZ, rotY, rotX] = rotate_cyls(angles);
        rotMat = rotZ* rotY* rotX;

        % Apply rotation to each point in the cylinder
        for i = 1:size(cylVertices, 2)
            for j = 1:size(cylVertices, 3)
                point = squeeze(cylVertices(:, i, j)) - center';
                rotatedPoint = rotMat * point;
                cylVertices(:, i, j) = rotatedPoint + center';
            end
        end

        % Update the cylinder with rotated points
        cyls2{electrodeNum} = cylVertices;

        % Redraw cylinder

        surf(squeeze(cylVertices(1, :, :)), squeeze(cylVertices(2, :, :)), squeeze(cylVertices(3, :, :)));
        hold on;

        c = c+1;
        rot_angles(c,:) = angles;
    end
end

% change rotation of ground electrode shank
cylVertices = cyls{32};

center = elec_coord_mm(32,:);
angles = mean(rot_angles,1, 'omitnan');
[rotZ, rotY, rotX] = rotate_cyls(angles);
rotMat = rotZ* rotY* rotX;
for i = 1:size(cylVertices, 2)
    for j = 1:size(cylVertices, 3)
        point = squeeze(cylVertices(:, i, j)) - center';
        rotatedPoint = rotMat * point;
        cylVertices(:, i, j) = rotatedPoint + center';
    end
end
cyls2{32} = cylVertices;
surf(squeeze(cylVertices(1, :, :)), squeeze(cylVertices(2, :, :)), squeeze(cylVertices(3, :, :)));

% Connecting the electrode shanks.

% scalar = 8;
%
% for l = 1:length(rod)
%
%     shank = elec_coord_mm(rod{l},:)';
%
%     if size(shank,2) == 1
%         p1 = shank(:,end) + [0, 0, 1.2]'; p2 = p1 - [0, 0, 2.4]';
%         scalar = 30;
%     else
%         p1 = shank(:,end); p2 = shank(:,end-1);
%     end
%
%     distance = [(p1(1) - p2(1)), (p1(2) - p2(2)), (p1(3) - p2(3))];
%
%     vec = (distance * scalar) + p1';
%
%     shank = [shank, vec'];
%     f = cscvn(shank(:,1:end));
%     fnplt(f, 10) %,'linewidth', 2);
%     hold on
%
%     shank_spline{l} = f;
% end


% Connecting the electrode shanks.

scalar = 8;

for l = 1:length(rod)

    num_elecs = size(rod{l},1);

    node_splines = [];

    % get the coordinates of all cylinders along each outer node of the
    % cylinger to connect these
    for s = 1:length(rod{l})
        if num_elecs == 1
            node_splines(1,:,:) = squeeze(cyls2{rod{l}(s)}(:,5,:));
            node_splines(2,:,:) = squeeze(cyls2{rod{l}(s)}(:,11,:));

        else
            node_splines(s,:,:) = squeeze(cyls2{rod{l}(s)}(:,6,:));
        end

    end

    for p = 1:size(node_splines,3)
        shank = node_splines(:,:,p)';

        % if size(shank,2) == 1
        %     p1 = shank(:,end) + [0, 0, 1.2]'; p2 = p1 - [0, 0, 2.4]';
        %     scalar = 30;
        %
        % else
        p1 = shank(:,end); p2 = shank(:,end-1);
        % end

        distance = [(p1(1) - p2(1)), (p1(2) - p2(2)), (p1(3) - p2(3))];

        if num_elecs == 1
            scalar = 35;
        end

        vec = (distance * scalar) + p1';

        % if size(shank,2) == 1
        %     angles = mean(rot_angles,1, 'omitnan');
        %     [rotZ, rotY, rotX] = rotate_cyls(angles);
        %     rotMat = rotZ* rotY* rotX;
        %
        %     rotatedPoint = rotMat * vec';
        %     vec = rotatedPoint';
        % end

        shank = [shank, vec'];
        f = cscvn(shank(:,1:end));
        fnplt(f) %,'linewidth', 2);
        hold on

        shank_spline{l,p} = f;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make a binary mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to pixel space

transformation_mat_pixel = [[pixel_size(1); 0; 0], [0; pixel_size(2); 0], [0; 0; pixel_size(3)]];
t2 = inv(transformation_mat_pixel);
pix_elec_coord = elec_coord_mm * t2;

% new cylinders in pixel space


for shankIdx = 1:length(rod) % Iterate over each shank
    for cylIdx = 1:length(rod{shankIdx}) % Iterate over each cylinder in the shank
        electrodeNum = rod{shankIdx}(cylIdx);
        cylVertices = cyls2{electrodeNum}; % 3 x 12 x 50 matrix

        for node = 1:size(cylVertices,3)
            node_coord = squeeze(cylVertices(:,:,node));
            j = node_coord' * t2;
            pix_cyl(:,:,node) = j';
        end

        % Update the cylinder with rotated points
        cyls_pixel{electrodeNum} = pix_cyl;

        % Redraw cylinder

        surf(squeeze(pix_cyl(1, :, :)), squeeze(pix_cyl(2, :, :)), squeeze(pix_cyl(3, :, :)));
        hold on;

    end
end


% New splines in pixel space
%
% Connecting the electrode shanks.
%
% scalar = 16;
% for l = 1:length(rod)
%
%     shank = elec_coord_mm(rod{l},:)';
%
%     if size(shank,2) == 1
%         p1 = shank(:,end) + [0, 0, 1.2]'; p2 = p1 - [0, 0, 2.4]';
%         scalar = 60;
%     else
%         p1 = shank(:,end); p2 = shank(:,end-1);
%     end
%
%     distance = [(p1(1) - p2(1)), (p1(2) - p2(2)), (p1(3) - p2(3))];
%
%     vec = (distance * scalar) + p1';
%
%
%     shank = [shank, vec'];
%
%     shank2 = shank' * t2;
%     shank2 = shank2';
%     f = cscvn(shank2(:,1:end));
%     fnplt(f, 10) %,'linewidth', 2);
%     hold on
%
%
%     pixel_spline{l} = f;
% end
%


scalar = 16;

for l = 1:length(rod)

    node_splines = [];

    num_elecs = size(rod{l},1);

    % get the coordinates of all cylinders along each outer node of the
    % cylinger to connect these
    for s = 1:length(rod{l})

        if num_elecs == 1
            node_splines(1,:,:) = squeeze(cyls2{rod{l}(s)}(:,5,:));
            node_splines(2,:,:) = squeeze(cyls2{rod{l}(s)}(:,11,:));
        else
            node_splines(s,:,:) = squeeze(cyls2{rod{l}(s)}(:,6,:));
        end
    end

    for p = 1:size(node_splines,3)
        shank = node_splines(:,:,p)';


        p1 = shank(:,end); p2 = shank(:,end-1);

        distance = [(p1(1) - p2(1)), (p1(2) - p2(2)), (p1(3) - p2(3))];

        if num_elecs == 1
            scalar = 35;
        end

        vec = (distance * scalar) + p1';

        shank = [shank, vec'];
        shank2 = shank' * t2;
        shank2 = shank2';
        f = cscvn(shank2(:,1:end));
        fnplt(f) %,'linewidth', 2);
        hold on

        pixel_spline{l, p} = f;
    end
end

%% Convert pixels to binary mask

% clearvars -except cyls_pixel rod pix_cyl pixel_spline pixel_dimensions

mask_volume = zeros(pixel_dimensions);

% Mark the voxels in the mask corresponding to cylinders and splines as 1
for shankIdx = 1:length(rod)
    % For cylinders
    for cylIdx = 1:length(rod{shankIdx})
        electrodeNum = rod{shankIdx}(cylIdx);
        pix_cyl = cyls_pixel{electrodeNum};

        % Convert pixel coordinates to voxel indices
        voxel_indices = round(pix_cyl);

        % Mark the voxels corresponding to the cylinder as 1 in the mask
        mask_volume(voxel_indices(1, :, :), voxel_indices(2, :, :), voxel_indices(3, :, :)) = 1;

        disp(size(mask_volume))
    end

    % For splines

    for p = 1:size(node_splines,3)
        shank_spline = pixel_spline{shankIdx, p};
        t = linspace(shank_spline.breaks(1), shank_spline.breaks(end), 1000);
        spline_points = fnval(shank_spline, t);

        % Convert pixel coordinates to voxel indices
        voxel_indices = round(spline_points);
        idx = find(voxel_indices > pixel_dimensions(1));
        [row, col] = ind2sub(size(voxel_indices), idx);

        if ~isempty(row)
            voxel_indices(:, col) = [];
        end

        % Mark the voxels corresponding to the spline as 1 in the mask
        for i = 1:size(voxel_indices, 2)
            mask_volume(voxel_indices(1, i), voxel_indices(2, i), voxel_indices(3, i)) = 1;
        end

    end

    disp(['shank number ', num2str(shankIdx), ' mask_vol_size = ', num2str(size(mask_volume))])
end


%% convert to binary as its much smaller! 

mask_volume = mask_volume > 0.1; 

%%

save(['volume_data_rescale_', num2str(512 *res_factor), '.mat'], 'mask_volume', '-v7.3')
