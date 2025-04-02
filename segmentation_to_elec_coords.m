
data = load('transformed_elec_coord.mat', 'scirunnrrd'); 
data = data.scirunnrrd; 

% Check the structure of 'axis_table'
disp('Structure of axis table:');
axis_table = data.axis; 
data_matrix = data.data; 

% Extract voxel size and min/max range from the axis table
voxel_size_x = axis_table(1).spacing;  % Voxel size in x
voxel_size_y = axis_table(2).spacing;  % Voxel size in y
voxel_size_z = axis_table(3).spacing;  % Voxel size in z

min_x = axis_table(1).min;  % Minimum x-coordinate in mm
min_y = axis_table(2).min;  % Minimum y-coordinate in mm
min_z = axis_table(3).min;  % Minimum z-coordinate in mm

% Step 1: Find the voxel indices where electrodes are (i.e., where data_matrix == 1)
[elec_voxels_x, elec_voxels_y, elec_voxels_z] = ind2sub(size(data_matrix), find(data_matrix == 1));

% Step 2: Convert voxel indices to real-world coordinates in mm
real_coords_x = min_x + (elec_voxels_x - 1) * voxel_size_x;
real_coords_y = min_y + (elec_voxels_y - 1) * voxel_size_y;
real_coords_z = min_z + (elec_voxels_z - 1) * voxel_size_z;

% Combine the real-world coordinates into a single matrix
electrode_coords_mm = [real_coords_x'; real_coords_y'; real_coords_z'];

% Step 3: Cluster the coordinates into 64 electrodes
% (This step assumes the electrodes are grouped relatively close together)

% If you have specific clustering criteria (e.g., spatial distance), you could apply a clustering algorithm like k-means.
% Here we will use k-means clustering with 64 clusters.
% num_electrodes = 49;
% [cluster_idx, w] = kmeans(electrode_coords_mm', num_electrodes);

% Parameters for DBSCAN
eps = 1;         % Maximum distance (in mm) between points to be considered as neighbors
minPts = 15;     % Minimum number of points to form a cluster (adjust based on your data)

% Perform DBSCAN clustering
cluster_idx = dbscan(electrode_coords_mm', eps, minPts);
num_electrodes = sum(unique(cluster_idx) > 0);

% Step 4: Find the coordinates of each electrode by taking the mean of the clustered points
electrode_centroids = zeros(num_electrodes, 3);  % 64 electrodes, each with x, y, z coordinates

for i = 1:num_electrodes
    electrode_centroids(i, :) = mean(electrode_coords_mm(:,cluster_idx == i), 2);

    scatter3(electrode_coords_mm(1,cluster_idx == i), electrode_coords_mm(2,cluster_idx == i),electrode_coords_mm(3,cluster_idx == i))
    hold on 
end



% Optionally, save the electrode coordinates to a CSV file
csvwrite('electrode_coordinates_mm.csv', electrode_centroids);


