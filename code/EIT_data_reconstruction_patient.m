%% Creat EIT signal recon image VTK files 

clc 
clear 
close all 

%% Preamble 
best_mesh = 13; 

inpath = '/home/eit/Mesh_convergence/meshes_to_test_pert'; 

load([inpath, '/names_to_test.mat']); 

load([name_test{best_mesh},'img.mat']); 

load('/home/yuval/PhD/Spike_detection/E05/EIT003/EIT_spikes_new.mat')

run('/home/eit/Packages_and_tools/eidors-v3.10-ng/eidors/startup.m')
clc

addpath('/home/eit/Packages_and_tools/Random_functions')
addpath('/home/eit/Packages_and_tools/Reconstruction-master/')

%% Calculate forward solution and jacobian 

v1 = fwd_solve(img); 

J = calc_jacobian(img); 

%% Create VTKs for the EIT signal for cluster 1

 
[Mesh_hex,J_hex] = convert_fine2coarse(img.fwd_model.elems,img.fwd_model.nodes,J,4 * 1e-3); % make element size 0.5mm

cd(inpath)

if ~exist('recon_images','dir')
    mkdir('recon_images')
end 

cd('recon_images');

cluster1 = squeeze(clust_EIT_spikes_raw(1,:,:));
figure; 
plot(cluster1)
title('All Data Cluter') 


rm_chan = std(cluster1, [], 1); 
cluster1_clean = cluster1; 
cluster1_clean(:, rm_chan > 3) = []; 

figure; 
plot(cluster1_clean)
title('Noisey chans removed Data Cluter') 

cluster1_bs = (cluster1_clean - mean(cluster1_clean(1:50,:))); 

figure; 
plot(cluster1_bs)
title('Baseline Normalised Data Cluter') 


d = v1.meas; 
d(v1.meas < 0) = -1; 
d(v1.meas > 0) = 1; 

d(rm_chan > 3) = []; 

flip = abs(cluster1_clean); 
flip = flip .* d'; 

J_hex_clean = J_hex(rm_chan <3,:); 

[Sigma,X] = eit_recon_tik0(flip,J_hex_clean,logspace(-25,2,250));

timpnts = [25:40:size(Sigma,2), 200]; 

for t =1:length(timpnts)% 60:70 %size(Sigma,2)
    
    ss = Sigma(:,[timpnts(t) - 10: timpnts(t) + 10]); 
    ss = mean(ss,2); 
    
    imgr = img;
    
    for i = 1:length(Mesh_hex.cells)
        imgr.elem_data(Mesh_hex.cells{i}) = ss(i);
    end
    
    
    writeVTKcell_orig(['EIT_recon_signal_t',num2str(timpnts(t))],imgr.fwd_model.elems,imgr.fwd_model.nodes,imgr.elem_data)
    
    
    
end
save('Recon_EIT_signal', 'Sigma');
    


%% Actual Spikes location 


true_spike = [0.146710206513164	0.0890890640195173	0.133260819149390]; 
c = find_element_centres(img.fwd_model);
 
idx = find_perturbation_indices(c,1,0.005,true_spike);

imgi = img; 

imgi.elem_data(idx{1}) = img.elem_data(idx{1}) .* 1; %.010;     

v2 = fwd_solve(imgi);

dv = v2.meas - v1.meas;
        
[Sigma4,X2] = eit_recon_tik0(dv',J_hex,logspace(-25,2,250));

imgr = imgh;

for i = 1:length(Mesh_hex.cells)
    imgr.elem_data(Mesh_hex.cells{i}) = Sigma2(i);
end


 
 
