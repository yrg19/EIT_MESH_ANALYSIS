%% add surface electrodes

addpath(genpath('/home/yuval/Upgrade_report'))
addpath(genpath('/home/eit/Packages_and_tools/Random_functions'))
run('/home/eit/Packages_and_tools/eidors-v3.10-ng/eidors/startup.m')

%%
Nasion = extract_fcsv('cropped_nasion.fcsv'); 
Left = extract_fcsv('cropped_left.fcsv'); 
Right = extract_fcsv('cropped_right.fcsv'); 

%% 

cd('/home/yuval/MATLAB Add-Ons/Collections/FieldTrip')
ft_defaults
clc 

% Define the 10-20 system electrode positions
cfg = [];
cfg.elec = ft_read_sens('standard_1020.elc'); % Standard 10-20 positions file

% Create a structure for electrodes
elec = ft_convert_units(cfg.elec, 'mm');  

electrode_labels = strtrim(string(elec.label));

% Define the standard 32 gray electrodes from the image
standard_32_labels = { ...
    'Fp1', 'Fp2', ...
    'AF3', 'AF4', ...
    'F7', 'F3', 'Fz', 'F4', 'F8', ...
    'FC1', 'FC2', 'FC6', 'FC5', ...
    'T7', 'C3', 'Cz', 'C4', 'T8', ...
    'CP5', 'CP6', 'CP1', 'CP2', ...
    'P7', 'P3', 'Pz', 'P4', 'P8', ...
    'PO3', 'PO4', ...
    'O1', 'Oz', 'O2'};


% Extract only the standard 32 electrodes
standard_32 = [];
for i = 1:length(standard_32_labels)
    idx = find(strcmp(electrode_labels, standard_32_labels{i}));
    if ~isempty(idx)
        standard_32.label{i} = elec.label{idx};
        standard_32.elecpos(i, :) = elec.elecpos(idx, :);
    end
end


% Define fiducial points for electrodes and scalp
elec_fiducials = [elec.elecpos(strcmp(elec.label, 'Nz'), :);
                  elec.elecpos(strcmp(elec.label, 'LPA'), :);
                  elec.elecpos(strcmp(elec.label, 'RPA'), :)];

scalp_fiducials = [Nasion; Left; Right];

% Estimate transformation using Procrustes analysis
[d, Z, transform] = procrustes(scalp_fiducials, elec_fiducials, 'Reflection',false);

% Apply transformation to all electrode positions
aligned_elecpos = transform.b * standard_32.elecpos * transform.T + transform.c(1,:);

scatter3_lazy(aligned_elecpos, 'r')
hold on
scatter3_lazy(standard_32.elecpos, 'g')
scatter3_lazy(elec_fiducials, 'b')
scatter3_lazy(scalp_fiducials, 'k')



%% Create model

fmdl = mk_common_model('a2C');
fmdl = fmdl.fwd_model;

%% Generate mesh 

inpath = ('/home/eit/Mesher/output');

name = 'Patient1_test_'; 

Tetra = readmatrix([inpath, '/', name,'tetra.csv']);
                                     
Nodes = readmatrix([inpath, '/', name,'vertices.csv']);


%% Align scalp electrodes 

fmdl.elems = Tetra(Tetra(:,5) == 1, 1:4); %

fmdl.nodes = Nodes;

fmdl.boundary = dubs3_2(fmdl.elems);

fmdl.electrode = struct('nodes',[],'z_contact',[]);

sign = ones(sum(Tetra(:,5) == 1),1) * 0.44; % ones(length(Tetra), 1) .* 0.44; %


%% Align Electrodes to scalp 

boundary_nodes = fmdl.nodes(unique(fmdl.boundary(:)), :); 
aligned_elecpos_m = aligned_elecpos ./ 1000; 

% Loop through each electrode and find the nearest node on the scalp mesh
for i = 1:size(aligned_elecpos_m, 1)
    % Compute distance from the current electrode to all nodes in the mesh
    distances = sqrt(sum((boundary_nodes - aligned_elecpos_m(i, :)).^2, 2));
    
    % Find the closest node
    [~, min_idx] = min(distances);
    
    % Update the electrode position to the closest node
    corrected_elecpos(i, :) = boundary_nodes(min_idx, :);
end

figure
scatter3_lazy(corrected_elecpos)

%% Create Model of All Tissues and add Scalp + Depth electrodes 

Tetra_orig = Tetra;
Tetra(Tetra(:,5) == 7,:) = [];

fmdl.elems = Tetra(:,1:4);

fmdl.nodes = Nodes;

fmdl.boundary = dubs3_2(fmdl.elems);

fmdl.electrode = struct('nodes',[],'z_contact',[]);

sig = zeros(length(Tetra),1);
%sig(Tetra(:,5) == 7) = 0.00001; % electrodes
sig(Tetra(:,5) == 1) = 0.44; % scalp
sig(Tetra(:,5) == 2) = 1.79; % csf
sig(Tetra(:,5) == 3) = 0.0001; % air
sig(Tetra(:,5) == 4) = 0.018; % skull
sig(Tetra(:,5) == 5) = 0.15; % white matter
sig(Tetra(:,5) == 6) = 0.3; % grey matter

sign = sig; 

%% Show mesh surface with scalp electrodes to ensure correct positions 

depth_elecs = readmatrix('electrode_coordinates_mm.txt'); 
depth_elecs(:,3) = depth_elecs(:,3) - 30.223092; % adjust for cropping (calculated by multiplying pixle size by number of pixles cropped)
depth_elecs_m = depth_elecs ./ 1000; 

img = mk_image(fmdl,sign');

img = add_electrode(img,[corrected_elecpos; depth_elecs_m],0.0015,0.2); 

figure
show_fem(img.fwd_model)

%% Save Electrode Coordinates 

all_elecs = [corrected_elecpos; depth_elecs_m]; 
all_elecs_mm = all_elecs .* 1000; 

transformation_mat_pixel = [[ps(1); 0; 0], [0; ps(2); 0], [0; 0; ps(3)]];
t2 = inv(transformation_mat_pixel);
pix_coord = all_elecs_mm * t2;

figure; 
subplot(1,2,1)
scatter3_lazy(all_elecs)

subplot(1,2,2)
scatter3_lazy(pix_coord)

%% Save Pixel COORDS

cd('/home/yuval/Upgrade_report/Mesh_optimisation/Segmentations/')

% Open a text file for writing
fileID = fopen('all_coordinates_pix.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(pix_coord, 1)
    fprintf(fileID, '%.3f,%.3f,%.3f\n', pix_coord(i, 1), pix_coord(i, 2), pix_coord(i, 3));
end

% Close the file
fclose(fileID);

%% Save mm COORDS

% Open a text file for writing
fileID = fopen('all_coordinates_mm.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(all_elecs_mm, 1)
    fprintf(fileID, '%.3f,%.3f,%.3f\n', all_elecs_mm(i, 1), all_elecs_mm(i, 2), all_elecs_mm(i, 3));
end

% Close the file
fclose(fileID);

%% Save m COORDS

% Open a text file for writing
fileID = fopen('all_coordinates_m.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(all_elecs, 1)
    fprintf(fileID, '%.3f,%.3f,%.3f\n', all_elecs(i, 1), all_elecs(i, 2), all_elecs(i, 3));
end

% Close the file
fclose(fileID);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

inp_file_path = ('/home/yuval/Upgrade_report/Mesh_optimisation/Segmentations');

load([inp_file_path,'/transformed_elec_coord.mat'])

pd = [scirunnrrd.axis(1).size,scirunnrrd.axis(2).size, scirunnrrd.axis(3).size];

ps = [scirunnrrd.axis(1).spacing,scirunnrrd.axis(2).spacing, scirunnrrd.axis(3).spacing];


%% 

depth_elecs_pix = readmatrix('electrode_coordinates_pix.txt'); 

depth_elecs = readmatrix('electrode_coordinates_mm.txt'); 
depth_elecs_m = depth_elecs ./ 1000; 

corrected_elecpos_mm = corrected_elecpos .* 1000; 

scatter3(depth_elecs(:,1), depth_elecs(:,2), depth_elecs(:,3))

hold on 

scatter3(corrected_elecpos_mm(:,1), corrected_elecpos_mm(:,2), corrected_elecpos_mm(:,3))

hold off 

transformation_mat_pixel = [[ps(1); 0; 0], [0; ps(2); 0], [0; 0; ps(3)]];
t2 = inv(transformation_mat_pixel);
pix_scalp_coord = corrected_elecpos_mm * t2;

transformation_mat_pixel = [[ps(1); 0; 0], [0; ps(2); 0], [0; 0; ps(3)]];
t2 = inv(transformation_mat_pixel);
pix_depth_coord = depth_elecs * t2;

figure

scatter3(pix_depth_coord(:,1), pix_depth_coord(:,2), pix_depth_coord(:,3))

hold on 

scatter3(pix_scalp_coord(:,1), pix_scalp_coord(:,2), pix_scalp_coord(:,3))

hold off 

figure

scatter3(pix_depth_coord(:,1), pix_depth_coord(:,2), pix_depth_coord(:,3))

hold on 

scatter3(depth_elecs_pix(:,1), depth_elecs_pix(:,2), depth_elecs_pix(:,3))

%% 

all_coords = [pix_depth_coord; pix_scalp_coord]; 

% Open a text file for writing
fileID = fopen('all_coordinates_pix.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(all_coords, 1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\n', all_coords(i, 1), all_coords(i, 2), all_coords(i, 3));
end

% Close the file
fclose(fileID);

%% 

all_coords = [depth_elecs; aligned_elecpos]; 

% Open a text file for writing
fileID = fopen('all_coordinates_mm.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(all_coords, 1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\n', all_coords(i, 1), all_coords(i, 2), all_coords(i, 3));
end

% Close the file
fclose(fileID);

%% 

all_coords = [depth_elecs; aligned_elecpos]; 

all_coords = all_coords ./ 1000; 

% Open a text file for writing
fileID = fopen('all_coordinates_m.txt', 'w'); % 'w' means write mode

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file for writing.');
end

% Write the coordinates to the text file, one row per line
for i = 1:size(all_coords, 1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\n', all_coords(i, 1), all_coords(i, 2), all_coords(i, 3));
end

% Close the file
fclose(fileID);