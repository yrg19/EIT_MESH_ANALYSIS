%% Scalp electrodes effect 


clear 
close all 
clc

%% Preamble 

cd('/home/yuval/Packages_and_tools/eidors-v3.10-ng/eidors')
startup
clc

addpath('/home/yuval/Packages_and_tools/Random_functions')


inpath = ('/home/yuval/Documents/Upgrade_report/Mesh_convergence/Meshes');

cd(inpath)

if ~exist('scalp_elecs', 'dir')
    mkdir('scalp_elecs')
end 

outpath = ('/home/yuval/Documents/Upgrade_report/Mesh_convergence/Meshes/scalp_elecs');

%% Load mesh 

cd(inpath)

mesh_name = 'Patient1_mesh_10_coarse_0.5_fine_2_rad_5_ratio_15_angle_1_dist_';

load([mesh_name,'img'])

n_elec = length(img.fwd_model.electrode); 

reorder = [img.fwd_model.electrode(33:n_elec), img.fwd_model.electrode(1:32)]; 

% change it so that depth electrodes are first 
img.fwd_model.electrode = reorder; 

for i = 1:49 %length(img.fwd_model.electrode)

    scatter3_lazy(img.fwd_model.nodes(img.fwd_model.electrode(i).nodes,:))

    hold on 

    p = img.fwd_model.nodes(img.fwd_model.electrode(i).nodes(1),:); 

    text(p(1), p(2), p(3) + 0.002, num2str(i))

end 

reference = 19; 

%% Make protocol without scalp elecs 


depth = [1:18, 20:49];

% Randomly shuffle the numbers
shuffled_numbers = depth(randperm(length(depth)));

% Select the first 40 numbers (20 pairs of unique numbers)
selected_numbers = shuffled_numbers(1:48);

% Reshape into 20 rows and 2 columns
random_pairs_no_repeats = reshape(selected_numbers', [24, 2]);

% Make stimulation protocol 
stim_depth = mk_stim_kai(random_pairs_no_repeats, length(depth) +1, 1e-4, 1,reference); 


% Make stimultion protocol with recording at scalp electrodes 
stim_scalp = mk_stim_kai(random_pairs_no_repeats, n_elec, 1e-4,1,reference); 


%% Find forward solution for both stimulation protocols 

img_orig = img; 

img.fwd_model.stimulation = stim_depth;
img.fwd_model.electrode = img.fwd_model.electrode(1:49); 
                                      
depth_forward = fwd_solve(img);
                        
img = img_orig; 
img.fwd_model.stimulation = stim_scalp;
                                      
scalp_forward = fwd_solve(img);

save([outpath,'/scalp_pre'], "scalp_forward", "depth_forward", "img", "stim_scalp", 'stim_depth', 'random_pairs_no_repeats')

%% %% Dedicing perturbation indices

elecs = readmatrix('/home/yuval/Documents/Upgrade_report/Input_files/all_coordinates_m.txt');

elecs2 = elecs(33:end,:);

% Parameters for DBSCAN
eps = 10e-3;         % Maximum distance (in mm) between points to be considered as neighbors
minPts = 4;     % Minimum number of points to form a cluster (adjust based on your data)

% Perform DBSCAN clustering
cluster_idx = dbscan(elecs2, eps, minPts);
num_electrodes = sum(unique(cluster_idx) > 0);

% Step 4: Find the coordinates of each electrode by taking the mean of the clustered points

electrode_centroids = zeros(num_electrodes, 3);  % 64 electrodes, each with x, y, z coordinates

colors = rand(num_electrodes, 3); % Nx3 matrix with values between 0 and 1

for i = 1:num_electrodes
    electrode_centroids(i, :) = mean(elecs2(cluster_idx== i, :), 1);

    scatter3_lazy(elecs2(cluster_idx == i, :), [colors(i,1), colors(i,2), colors(i,3)]) % , E05_elecs(2,cluster_idx == i),electrode_coords_mm(3,cluster_idx == i))
    hold on 
end

electrode_centroids = electrode_centroids .* 1000; 

%% Get random gitters to generate perturbations near electrodes 

% Create a 9x3 matrix (example)
A = rand(num_electrodes,3);  % Random values for demonstration

% Generate random values between 0.5 and 1
random_values = 3 + (1-0.5) * rand(num_electrodes,3);

% Generate random signs (-1 or 1)
random_signs = randi([0,1], num_electrodes,3) * 2 - 1; % Converts 0s to -1s and 1s remain 1s

% Apply random addition/subtraction
A_modified = A + (random_values .* -1); % .* random_signs);


%% Choose perturbation locations 


ele_near = (electrode_centroids + A_modified) ./1000; % add gitters and convert to meters 
ri = randperm(num_electrodes, 5); 

ele_med = (electrode_centroids + (A_modified + [-10 17 -15]))./1000; % add bigger gitters and convert to meters 
ele_far = abs(electrode_centroids + (A_modified + [30 20 20]))./1000;
ele_far(1:3,:) = ele_far(1:3,:) +  ([30 0 0]./1000);

pert_test = {ele_near(1,:), ele_near(2,:), ele_near(3,:), ele_med(1,:), ele_med(2,:)...
    ele_med(3,:), ele_med(6,:), ele_med(7,:), ele_far(1,:), ele_far(2,:), ele_far(4,:), ele_far(5,:), ele_far(7,:)};

pert_name = {'near1', 'near2', 'near3', 'med1', 'med2', 'med3','med4','med5', 'far1', 'far2', 'far3', 'far4', 'far5'};


%% Plot perturbations 

c = find_element_centres(img.fwd_model);

ch = [8, 12];
for p = 1:length(ch)

    figure 

    idx = find_perturbation_indices(c,1,0.005,pert_test{ch(p)});

    imgi=img;

    imgi.elem_data(idx{1}) = 10;

    show_fem(imgi)

    title(pert_name{ch(p)})

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end 

%% Save perturbations when happy 

cd(outpath)
save('Perturbations', 'pert_name', 'pert_test', 'A_modified', 'ele_near', 'electrode_centroids', 'ele_far', 'ele_med', 'ri'); 

%% Simulation reconstruction error

cd(outpath) 

load('Perturbations.mat')

img = img_orig; 

img_scalp = img; 
img_scalp.fwd_model.stimulation = stim_scalp;

img_depth = img; 
img_depth.fwd_model.electrode = img.fwd_model.electrode(1:49); 
img_depth.fwd_model.stimulation = stim_depth;

                        
calc_these = {img_depth, img_scalp}; 
fwds = {depth_forward, scalp_forward}; 
cond_name = {'Depth_Only', 'Scalp+Depth'}; 

for n = length(calc_these)

    disp(['Began Calculating ', cond_name{n}, ' Condition WOO'])

    img = calc_these{n}; 

    c = find_element_centres(img.fwd_model);

    J = calc_jacobian(img);

    imgh = img;
       
    v1 = fwds{n};

    [Mesh_hex,J_hex] = convert_fine2coarse(imgh.fwd_model.elems,imgh.fwd_model.nodes,J,4 * 1e-3); % make element size 0.5mm

    for p = 1:length(pert_name)
 
        
        idx = find_perturbation_indices(c,1,0.005,pert_test{p});
        
        imgi=img;    
            
        imgi.elem_data(idx{1}) = imgi.elem_data(idx{1}) .* 1.01;     

        v2 = fwd_solve(imgi);

        dv = v2.meas - v1.meas;
        
        noise_size=[200, length(dv)];
        Noise = 1*1e-6*randn(noise_size);

        [Sigma,X] = eit_recon_tik0(dv',J_hex,logspace(-25,2,250), Noise);

        imgr = imgh;

        for i = 1:length(Mesh_hex.cells)
            imgr.elem_data(Mesh_hex.cells{i}) = Sigma(i);
        end


        %near_elec.J = J;
        pert.Sigma{p}= Sigma;
        pert.X{p} = X;
        pert.dv{p} = dv;
        pert.imgi{p} = imgi; % unregularised pert
        pert.imgr{p} = imgr; % hexagonised, regularised image
        pert.mesh_hex{1} = Mesh_hex;
        pert.J_hex{1} = J_hex;
        pert.fwd_sol = v1; 
 

        writeVTKcell_orig([cond_name{n},'_','recon_mesh_', '_',pert_name{p}],imgr.fwd_model.elems,imgr.fwd_model.nodes,imgr.elem_data)

        writeVTKcell_orig([cond_name{n},'_','simulated_mesh_', '_',pert_name{p}],imgi.fwd_model.elems,imgi.fwd_model.nodes,imgi.elem_data)

        disp(['******************* Finished Calculating ', cond_name{n}, ' Perturbation ', num2str(p), ' ************************************'])
    
    end

    save([cond_name{n},'recon.mat'], "pert", '-v7.3')

    close all 

end

save('Scalp+Depthrecon.mat' , 'img_scalp' , 'scalp_forward', "-append")

save('Depth_Onlyrecon.mat' , 'img_depth' ,"-append")


%% Calculate WSV for both conditions 


WSV_true = zeros(length(cond_name), length(pert_test)); 
WSV_recon = zeros(length(cond_name), length(pert_test)); 

cd(outpath)
for n = 1:length(cond_name) 

    load([cond_name{n}, 'recon.mat'])

    for p = 1:length(pert_test)

        true_com = pert_test{p}; % true pert centre of mass 
    
        Sigma = pert.Sigma{p};

        Sigma(isnan(Sigma)) = 0;

        [WSV_true(n,p), WSV_recon(n,p)] = calc_wsv(img, pert.imgi{p},true_com, Sigma, pert.mesh_hex{1});

    end

    disp(['******* Calclated WSV for ', cond_name{n}, ' ***********'])

end

save('WSV_D_vs_SD', 'WSV_recon', 'WSV_true', 'pert_test', 'pert_name')

%% Plot Results

cond_name_plot = {'Depth Only', 'Scalp + Depth'}; 

reg_WSV = WSV_recon./WSV_true;

depth_WSV = reg_WSV(1,:);
scalp_WSV = reg_WSV(2,:);


scatter(1:length(pert_test), depth_WSV, 150, 'filled')
hold on
scatter(1:length(pert_test), scalp_WSV, 150, 'filled')
hold off;


ylabel('Corrected WSV', 'FontSize',16)
%ylim([0 20])
xlabel('Pertubation', 'FontSize',16)
xticks(1:length(pert_name))
xticklabels(pert_name)
title('WSV by Perturbation', 'FontSize',20)
legend(cond_name_plot, 'FontSize',16)

%% Determine influence of every scalp electrode

% load('Scalp+Depthrecon.mat')
% load('Perturbations.mat')

n_elec = length(img_scalp.fwd_model.electrode);

scalp_elecs = 50 -3 :n_elec -3; % scalp electrodes removing injection and reference

Mesh_hex = pert.mesh_hex{1};
J_hex = pert.J_hex{1};

c = find_element_centres(img_scalp.fwd_model);


%for e = 1:n_elec

sin_elec = [];


v1 = [];
   
e_img = img_scalp;

%     for s = 1:length(e_img.fwd_model.stimulation) % for every injection pair remove measurement through scalp electrode
%
%         mes = full(e_img.fwd_model.stimulation(s).meas_pattern);
%
%         mes(scalp_elecs(e),:) = 0;
%
%         e_img.fwd_model.stimulation(s).meas_pattern = mes;
%
%     end
%
%
v1 = scalp_forward;


for p = 11:length(pert_name)

    v2 = []; idx = []; imgi = []; imgr = []; dv = []; Noise = []; Sigma = []; X = [];

     
    disp([' --------------------- STARTING PERTURBAATION ', num2str(p), ' ----------------------'])
    
    idx = find_perturbation_indices(c,1,0.005,pert_test{p});

    imgi=e_img;

    imgi.elem_data(idx{1}) = imgi.elem_data(idx{1}) .* 1.01;

    v2 = fwd_solve(imgi);

    for e = 1:length(scalp_elecs)

       

        id_rm = scalp_elecs(e): 78: size(J_hex,1); % remove all electrode measurements

        f1 = v1.meas; f2 = v2.meas;

        f1(id_rm) = [];
        f2(id_rm) = [];

        dv = f2 - f1;

        noise_size=[200, length(dv)];
        Noise = 1*1e-6*randn(noise_size);

        J_hex2 = J_hex;
        J_hex2(id_rm,:) = [];

        [Sigma,X] = eit_recon_tik0(dv',J_hex2,logspace(-25,2,250), Noise);

        imgr = e_img;

        for i = 1:length(Mesh_hex.cells)
            imgr.elem_data(Mesh_hex.cells{i}) = Sigma(i);
        end


        sigma{e}= Sigma;
        x{e} = X;
        Dv{e} = dv;
        Imgi{e} = imgi; % unregularised pert
        Imgr{e} = imgr; % hexagonised, regularised image


        disp(['******* Done Calculating Perturbation ', num2str(p), ' of ',...
            num2str(length(pert_name)), ', Electrode ', num2str(e), ' of ', num2str(length(scalp_elecs)), ' removed ************'])
    end

    sin_elec.Sigma = sigma;
    sin_elec.X = x;
    sin_elec.dv = Dv;
    sin_elec.imgi = Imgi;
    sin_elec.imgr = Imgr;
    sin_elec.fwd_sol = v1;

    save(['sin_elec_pert', num2str(p), '_recon'], 'sin_elec', '-v7.3')

end

%% Calculate the WSV scores per pert and elec combo 


WSV_true_sin_elec = zeros(length(scalp_elecs), length(pert_test)); 
WSV_recon_sin_elec = zeros(length(scalp_elecs), length(pert_test)); 

cd(outpath)


for p = 1:length(pert_name) 

    load(['sin_elec_pert', num2str(p), '_recon.mat'])

    for e = 1:length(scalp_elecs)

        true_com = pert_test{p}; % true pert centre of mass 
    
        Sigma = sin_elec.Sigma{e}; 
    
        Sigma(isnan(Sigma)) = 0; 
    
        [WSV_true_sin_elec(e,p), WSV_recon_sin_elec(e,p)] = calc_wsv(img_scalp, sin_elec.imgi{p},true_com, Sigma, pert.mesh_hex{1}); 

    end 

    disp(['******* Calclated WSV for Perturbation', num2str(p), ' ***********'])

end 

save('WSV_sin_elec', 'WSV_recon_sin_elec', 'WSV_true_sin_elec', 'pert_test', 'pert_name')

 %% Plot results 

 load('WSV_D_vs_SD')

 reg_WSV = WSV_recon./WSV_true;

 scalp_WSV = reg_WSV(2,:);

 elec_WSV = WSV_recon_sin_elec ./ WSV_true_sin_elec;
% 
% % Add transition to soft yellow & peach
% cmap(:,2) = linspace(0.9, 0.8, 256)'; % Green channel fades slightly
% cmap(:,1) = linspace(0.8, 1, 256)';   % Add warmth (peach)
% cmap(:,3) = linspace(1, 0.9, 256)';   % Soften blue into lavender
% 
% 
% distributionPlot(elec_WSV, 'xNames', pert_name, 'showMM', 5, 'colormap', cmap);
% hold on 
% 
% plot(scalp_WSV, 'k', 'LineWidth', 2)
% 
% ylabel('Normalised WSV Values')

%% Find best electrodes 

ss = sum(elec_WSV,2); 

[ss_sort, id_ss] = sort(ss); 

figure()
plot(ss_sort)
xlabel('Electrode', 'FontSize',16)
ylabel('Sum of normalised WSV', 'FontSize',16)
hold on 
yline(sum(scalp_WSV))
legend({'Electrode Sum', ' Whole Scalp Sum'}, 'FontSize',14)
title('Effect of Electrode Removal', 'FontSize',18)


keep_elecs = id_ss(25); % these are the 6 most important electrodes
keep_elecs_real = keep_elecs + 49; % these are their actual numbers from 1:81

figure() 
for i = 1:n_elec

    coords = img_scalp.fwd_model.nodes(img_scalp.fwd_model.electrode(i).nodes,:); 

    if ismember(i, keep_elecs_real)
        scatter3_lazy(coords, 'g')
    else
        scatter3_lazy(coords, 'k')
    end

    hold on 
end

title('Important Scalp Electrodes')

% %% Kept Electrodes 
% 
% n_elec = length(img_scalp.fwd_model.electrode); 
% 
% scalp_elecs = 50 -3 :n_elec -3; % scalp electrodes removing injection and reference 
% 
% Mesh_hex = pert.mesh_hex{1}; 
% J_hex = pert.J_hex{1}; 
% 
% c = find_element_centres(img_scalp.fwd_model);
% 
% keep_elecs_meas = keep_elecs_real -4; % scalp electrode numbers after accounting for removal of 
% % injecting electrodes and reference 
% 
% e_img = img_scalp;
% 
% for s = 1:length(e_img.fwd_model.stimulation) % for every injection pair remove measurement through scalp electrode
% 
%     for e = 1:length(scalp_elecs)
% 
%         if ismember(scalp_elecs(e), keep_elecs_meas)
%             continue;
%         else
%             mes = full(e_img.fwd_model.stimulation(s).meas_pattern);
% 
%             mes(scalp_elecs(e),:) = 0;
% 
%             e_img.fwd_model.stimulation(s).meas_pattern = mes;
%         end
% 
%     end 
% end 
% 
% v1 = fwd_solve(e_img);
% 
% 
% for p = 1:length(pert_name)
% 
%     v2 = []; idx = []; imgi = []; imgr = []; dv = []; Noise = []; Sigma = []; X = [];
% 
%     idx = find_perturbation_indices(c,1,0.005,pert_test{p});
% 
%     imgi=e_img;
% 
%     imgi.elem_data(idx{1}) = imgi.elem_data(idx{1}) .* 1.01;
% 
%     v2 = fwd_solve(imgi);
% 
%     dv = v2.meas - v1.meas;
% 
%     noise_size=[200, length(dv)];
%     Noise = 1*1e-6*randn(noise_size);
% 
%     [Sigma,X] = eit_recon_tik0(dv',J_hex,logspace(-25,2,250), Noise);
% 
%     imgr = e_img;
% 
%     for i = 1:length(Mesh_hex.cells)
%         imgr.elem_data(Mesh_hex.cells{i}) = Sigma(i);
%     end
% 
% 
%     sigma{p}= Sigma;
%     x{p} = X;
%     Dv{p} = dv;
%     Imgi{p} = imgi; % unregularised pert
%     Imgr{p} = imgr; % hexagonised, regularised image
% 
% 
%     disp(['******* Done Calculating Perturbation ', num2str(p), ' of ',...
%         num2str(length(pert_name)), ' ************'])
% end
% 
% sin_elec.Sigma = sigma;
% sin_elec.X = x;
% sin_elec.dv = Dv;
% sin_elec.imgi = Imgi;
% sin_elec.imgr = Imgr;
% sin_elec.fwd_sol = v1;
% 
% save('scalp_min_elec_recon', 'sin_elec', '-v7.3')
% 
% 
% %% Calculate WSV scores with minimal scalp protocol 
% 
% for p = 1:length(pert_test)
% 
%     true_com = pert_test{p}; % true pert centre of mass
% 
%     Sigma = sin_elec.Sigma{p};
% 
%     Sigma(isnan(Sigma)) = 0;
% 
%     [WSV_true_min_elec(p), WSV_recon_min_elec(p)] = calc_wsv(img_scalp, sin_elec.imgi{p},true_com, Sigma, pert.mesh_hex{1});
% 
% end
% 
% WSV_min_reg = WSV_recon_min_elec./ WSV_true_min_elec; 
% 
% %% Compare Min_scalp, Scalp + depth and depth for reconstruction 
% 
% 
% load('WSV_D_vs_SD')
% 
% reg_WSV = WSV_recon./WSV_true; 
% 
% scalp_WSV = reg_WSV(2,:); 
% 
% depth_WSV = reg_WSV(1,:); 
% 
% scatter(1:13,WSV_min_reg, 150,'filled')
% hold on 
% scatter(1:13,scalp_WSV, 150,'filled')
% hold on 
% scatter(1:13,depth_WSV, 150,'filled') 
% 
% legend({'Minimal Scalp + Depth', ' Scalp + Depth', ' Depth Only'}, 'FontSize', 16)
% xlabel('Perturbation', 'FontSize', 16)
% xticks(1:13)
% ylabel('Normalised WSV score', 'FontSize', 16)
% xticklabels(pert_name)
% title('Reconstruction Accuracy by Protocol', 'FontSize', 18)
% 

%% NEW COMPARISON ANALYSIS 

ss = sum(elec_WSV,2); 

[ss_sort, id_ss] = sort(ss); 

n_elec = length(img_scalp.fwd_model.electrode); 

scalp_elecs = 50 -3 :n_elec -3; % scalp electrodes removing injection and reference 

Mesh_hex = pert.mesh_hex{1}; 
J_hex = pert.J_hex{1}; 

c = find_element_centres(img_scalp.fwd_model);


v1 = scalp_forward;


for p = length(pert_name):-1:1

    v2 = []; idx = []; imgi = []; sim_elec = []; sigma = []; Dv = []; X = []; Imgr= []; 

    idx = find_perturbation_indices(c,1,0.005,pert_test{p});

    imgi=e_img;

    imgi.elem_data(idx{1}) = imgi.elem_data(idx{1}) .* 1.01;

    v2 = fwd_solve(imgi);

    id_rm = []; 

    for e = 1:length(scalp_elecs)

        imgr = []; dv = []; Noise = []; Sigma = []; X = [];

        rm_e = scalp_elecs(id_ss(e)); 

        id_rm = [id_rm, rm_e:78:size(J_hex,1)]; % all measurements using this channel 
 
        J_hex2 = J_hex; 

        J_hex2(id_rm,:) = []; 

        forward_1 = v1.meas; 
        forward_1(id_rm,:) =[]; 

        forward_2 = v2.meas; 
        forward_2(id_rm,:) = []; 

        dv = forward_2 - forward_1;

        noise_size=[200, length(dv)];
        Noise = 1*1e-6*randn(noise_size);

        [Sigma,X] = eit_recon_tik0(dv',J_hex2,logspace(-25,2,250), Noise);

        imgr = e_img;

        for i = 1:length(Mesh_hex.cells)
            imgr.elem_data(Mesh_hex.cells{i}) = Sigma(i);
        end


        sigma{e}= Sigma;
        x{e} = X;
        Dv{e} = dv;
        Imgr{e} = imgr; % hexagonised, regularised image


        disp(['******* Done Calculating Perturbation ', num2str(p), ' of ',...
            num2str(length(pert_name)), ' Electrodes Removed = ', num2str(e),' ************'])

    end 

sin_elec.Sigma = sigma;
sin_elec.X = x;
sin_elec.dv = Dv;
sin_elec.imgi = imgi;
sin_elec.imgr = Imgr;
sin_elec.fwd_sol = v1;

save(['cumulative_elec_removal_pert', num2str(p)], 'sin_elec', '-v7.3')

end 

%% Calculate WSV by cumulative electrode removal 

for p = 1:length(pert_test) 

    load(['cumulative_elec_removal_pert', num2str(p)])

    for e = 1:length(scalp_elecs)

        true_com = pert_test{p}; % true pert centre of mass 
    
        Sigma = sin_elec.Sigma{e}; 
    
        Sigma(isnan(Sigma)) = 0; 
    
        [WSV_true_cum_elec(e,p), WSV_recon_cum_elec(e,p)] = calc_wsv(img_scalp, sin_elec.imgi,true_com, Sigma, pert.mesh_hex{1}); 

    end 

    disp(['******* Calclated WSV for Perturbation', num2str(p), ' ***********'])

end 


save('WSV_cumulative', 'WSV_recon_cum_elec', "WSV_true_cum_elec")
%% Plot WSV cumulative scores 

cum_reg = WSV_recon_cum_elec ./ WSV_true_cum_elec; 

norm_reg = cum_reg - scalp_WSV; 


% scatter(1:length(scalp_elecs), norm_reg(:,[6, 8:12]), 100, 'filled')
% hold on 
% plot(1:length(scalp_elecs), norm_reg(:,[6, 8:12]))
% 
% 
% cum_reg = WSV_recon_cum_elec ./ WSV_true_cum_elec; 
% norm_reg = cum_reg - scalp_WSV; 


% Extract the required columns
data = norm_reg(:, [6, 8:12]);

% Get the number of datasets being plotted
numLines = size(data, 2);

% Generate colors
colors = lines(numLines); % MATLAB's "lines" colormap gives distinct colors

scatterHandles = gobjects(numLines, 1); % Preallocate for scatter handles

hold on 
for i = 1:numLines
    scatterHandles(i) = scatter(1:length(scalp_elecs), data(:, i), 100, colors(i, :), 'filled'); 
    plot(1:length(scalp_elecs), data(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
end
hold off


xlabel('Number of electrodes removed', 'FontSize',16)
ylabel('Normalised WSV - Scalp WSV', 'FontSize',16)
title('Cumulative electrode removal WSV relative to whole scalp protocol', 'FontSize',18)
% Extract the correct legend labels
legendLabels = pert_name([6, 8:12]); 
legend(scatterHandles,legendLabels, 'FontSize',16)


%% Make VTK files for reconstructions 


load('cumulative_elec_removal_pert8.mat')


writeVTKcell_orig('minimal_scalp_protocol_med5',sin_elec.imgr{1}.fwd_model.elems,sin_elec.imgr{16}.fwd_model.nodes,sin_elec.imgr{16}.elem_data)



%% Rererference Scalp 

cd(outpath) 

load('Perturbations.mat')
load('Scalp+Depthrecon.mat')


c = find_element_centres(img_scalp.fwd_model);

J_hex = pert.J_hex{1}; 

imgh = img_scalp;

v1 = scalp_forward;

J_reref = J_hex;

id_rr = [];
for s = 1:length(img_scalp.fwd_model.stimulation)
    id_rr = 51+ (78*(s-1)):78+ (78*(s-1));
    J_reref(id_rr,:) = J_hex(id_rr,:) - J_hex(51+ (78*(s-1)),:);
end


for p = 1:length(pert_name)

    idx = find_perturbation_indices(c,1,0.005,pert_test{p});

    imgi=img_scalp;

    imgi.elem_data(idx{1}) = imgi.elem_data(idx{1}) .* 1.01;

    v2 = fwd_solve(imgi);

    dv = v2.meas - v1.meas;

    noise_size=[200, length(dv)];
    Noise = 1*1e-6*randn(noise_size);

    [Sigma,X] = eit_recon_tik0(dv',J_reref,logspace(-25,2,250), Noise);

    imgr = imgh;

    for i = 1:length(Mesh_hex.cells)
        imgr.elem_data(Mesh_hex.cells{i}) = Sigma(i);
    end


    %near_elec.J = J;
    pert.Sigma{p}= Sigma;
    pert.X{p} = X;
    pert.dv{p} = dv;
    pert.imgi{p} = imgi; % unregularised pert
    pert.imgr{p} = imgr; % hexagonised, regularised image
    pert.mesh_hex{1} = Mesh_hex;
    pert.J_hex{1} = J_reref;
    pert.fwd_sol = v1;

    disp(['******************* Finished Calculating Perturbation ', num2str(p), ' ************************************'])

end

save('reref_scalp_recon.mat', "pert", '-v7.3')

