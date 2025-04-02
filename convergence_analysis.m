%% Analysise meshes for convergence 

clear 
clc
close all 

%% Preamble  

run('/home/eit/Packages_and_tools/eidors-v3.10-ng/eidors/startup.m')
addpath('/home/eit/Packages_and_tools/Random functions')
addpath(genpath('/home/eit/Packages_and_tools/Reconstruction-master/'))

coarse=[40 10 5];  
fine= [1 0.5 0.2]; 
radius=[4 2.5 2]; 
facet_ratio=[10 5 2]; 
facet_angle=[30 15 10]; 
facet_distance=[1 0.5 0.1]; 

param_titles = {'Coarse Cell Size'; 'Fine Cell Size'; 'Radius Around Electrodes';'Facet Ratio'; 'Facet Angle';...
    'Facet Distance'};

params = {coarse, fine, radius, facet_ratio, facet_angle, facet_distance};


param_nums = [40 10 5; 1 0.5 0.2; 4 2.5 2; 10 5 2; 30 15 10; 1 0.5 0.1]; 

load('/home/yuval/Upgrade_report/Mesh_optimisation/combi_forwards_all.mat')

clc

%% Understand why meshes have not been processed 

mesh_removed = No_mesh; 

for i = 1:length(mesh_removed)

    if isempty(mesh_removed{i})
        mesh_r_reason(i) = 0; 
    elseif contains(mesh_removed{i}, 'elec_problem')
        mesh_r_reason(i) = 1; 
    elseif contains(mesh_removed{i}, 'too_big')
        mesh_r_reason(i) = 2; 
    else
        mesh_r_reason(i) = 3; 
    end 

end 

mesh_r_reason(Mesh_size > 45e6) = 2; 
figure
hist(mesh_r_reason)
xticks(0:3)
xticklabels({'Good Mesh', 'Elec Problem', 'Too Big', 'No Mesh'})

%% Get the indices per mesh 


i = 0; 
for c = 1:length(coarse)
    for f = 1:length(fine)
        for r = 1:length(radius)
            for fr = 1:length(facet_ratio)
                for fa = 1:length(facet_angle)
                    for fd = 1:length(facet_distance)

                        i = i+ 1; 
                        all_ind_map(i,:) = [c, f, r, fr, fa, fd]; 
                        
                        combi_mesh_fwd(i,:) = combi_op_fwd(c, f, r, fr, fa, fd,:); 
                        
                        
                    end 
                end 
            end 
        end 
    end 
end 

% too_small = find(mesh_r_reason ~= 0); 
% 
% B = all_ind_map(too_small,:); 
% 
% linearIndices = sub2ind([3, 3, 3, 3, 3, 3], ...
%                         B(:, 1), B(:, 2), B(:, 3), B(:, 4), B(:, 5), B(:, 6));
% 
% 
% num_small = num_elem(linearIndices); 


h = sum(combi_mesh_fwd,2); % only keep meshes that ran fully 

B = all_ind_map(find(h ~= 0),:); % keep parameters for meshes that ran fully 
% 
% linearIndices = sub2ind([3, 3, 3, 3, 3, 3], ...
%                         B(:, 1), B(:, 2), B(:, 3), B(:, 4), B(:, 5), B(:, 6));


num_survived = Mesh_size(find(h ~= 0)); 

[sorted,id] = sort(num_survived, 'ascend'); % sort meshes by size 



%% Parameter Differences in measurements 


%p_na = ~any(isnan(all_ind_map')); 
%sur_ind_map = all_ind_map(p_na,:); 
sorted_ind_map = B(id,:); 

h2 = combi_mesh_fwd(find(h ~= 0 ),:);
G = sum(h2(id,:)'.^2);

figure;
names ={'Coarse','Fine','Radius','Ratio','Angle','Distance'};
for p = 1:6
    
    % Prepare data 
    
    [labels, id2] = sort(sorted_ind_map(:,p)); 
    values  = G(id2);
    
    
    % Get Unique Labels
    uniqueLabels = unique(labels);
    
    % Organize Data for Plotting
    groupedData = cell(length(uniqueLabels), 1);
    
    b = 0; 
    for i = length(uniqueLabels):-1:1
        
        b = b+1; 
        groupedData{b} = values(labels == uniqueLabels(i));
       
    end
    
    % Violin Plot using distributionPlot
        
    subplot(2,3,p)

   
    % Add transition to soft yellow & peach
    cmap(:,2) = linspace(0.9, 0.8, 256)'; % Green channel fades slightly
    cmap(:,1) = linspace(0.8, 1, 256)';   % Add warmth (peach)
    cmap(:,3) = linspace(1, 0.9, 256)';   % Soften blue into lavender

    
    distributionPlot(groupedData, 'xNames', uniqueLabels, 'showMM', 5, 'colormap',cmap);
    
    
    % Increase font size of axes labels, title, and ticks
    set(gca, 'FontSize', 14)  % Change 14 to your preferred size
    xlabel('Parameter Size', 'FontSize', 16) % Adjust axis label size
    ylabel('Forward Solution Sum', 'FontSize', 16)
    xticks(1:3)
    xticklabels(fliplr(param_nums(p,:)))
    title(names{p}, 'FontSize', 18) % Title size
    
    
   
end




%% Compare dvs by number of elements and by parameter

G2 = combi_mesh_fwd(find(h~= 0),:)'; 
meas_sort = G2(:,id); 



%% Mesh convergence by number of elements 


% Number of meshes per bin
factor = 11;

% Calculate the number of bins
numMeshes = numel(sorted);
numBins = ceil(numMeshes / factor);

% Initialize bin assignments
binAssignments = zeros(numMeshes, 1);

% Group the meshes into bins of 11
for i = 1:numBins
    startIdx = (i-1) * factor + 1; % Starting index for each bin
    endIdx = min(i * factor, numMeshes); % Ending index for each bin (handles the last bin)
    
    binAssignments(startIdx:endIdx) = i; % Assign the bin number to the appropriate meshes
end

% Now you can group the meshes by binAssignments
bins = cell(numBins, 1);
for i = 1:numBins
    bins{i} = sorted(binAssignments == i);
    bt(i) = bins{i}(1);
end


var = []; size_cat = []; 

for w = 1:length(unique(binAssignments))

    size_cat(w) = length(find(binAssignments == w)); 
    
end 

if sum(size_cat) ~= length(num_survived)
    error('Wrong Bin Sizes')
end 

min_size = min(size_cat); 

for b = 1:length(unique(binAssignments))

    mat = meas_sort(:, binAssignments == b); % get the measuremetns for all meshes from size category 

    %mat = mat .^2; 
    
    for n = 1:1000


        idx = randperm(size_cat(b),round(min_size/2)); 

        var(b,n) = mean(std(mat(:,idx), [], 2));

    end 

    var_ind_map{b} = sorted_ind_map(binAssignments == b, :); 
end 

figure; 


% Add transition to soft yellow & peach
cmap(:,2) = linspace(0.9, 0.8, 256)'; % Green channel fades slightly
cmap(:,1) = linspace(0.8, 1, 256)';   % Add warmth (peach)
cmap(:,3) = linspace(1, 0.9, 256)';   % Soften blue into lavender


distributionPlot(var', 'xNames', unique(binAssignments), 'showMM', 5, 'colormap',cmap);


% Increase font size of axes labels, title, and ticks
set(gca, 'FontSize', 14)  % Change 14 to your preferred size

%scatter(1:length(unique(bin)), var)
    
% hold on 
% hLine = plot(mean(var,2)); % Line plot, store the handle in hLine
% legend(hLine, 'Mean Variance') % Create a legend only for the line plot


xticks(1:length(unique(binAssignments)))
%xticklabels(round(bt ./ 1e6,1))
xlabel('Number of Elements, Millions',  'FontSize', 16)
ylabel('Mean Permuted Variance Across Measurements', 'FontSize', 16)
ylim([-0.5e-4 2.5e-4])
title('Mesh Convergence by Element Number', 'FontSize', 18)



%% Meshes to test 

name_test = []; var_type = []; 

ave= mean(var,2); 

test = [9]; 

test2 = find(binAssignments == test(:)); 

i = 0;

for b = 1:length(test)
    
    tmp = var_ind_map{test(b)};
    
    
    for t = 1:length(tmp)
        
        i = i+1;
        name_test{i} = ['Patient1_mesh_', num2str(param_nums(1, tmp(t,1))), ...
            '_coarse_', num2str(param_nums(2, tmp(t,2))), ...
            '_fine_', num2str(param_nums(3, tmp(t,3))), ...
            '_rad_', num2str(param_nums(4, tmp(t,4))), ...
            '_ratio_', num2str(param_nums(5, tmp(t,5))), ...
            '_angle_', num2str(param_nums(6, tmp(t,6))), ...
            '_dist_'];
        
        
        var_type(i,:) = [param_nums(1, tmp(t,1)), param_nums(2, tmp(t,2)), param_nums(3, tmp(t,3)),...
            param_nums(4, tmp(t,4)), param_nums(5, tmp(t,5)), param_nums(6, tmp(t,6))];
        
        sizes(i) = sorted(test2(i)); 
        
        
    end
end

%% Save meshes to test in a new folder to transfer

cd('/home/yuval/Upgrade_report/Mesh_optimisation')
if ~exist('meshes_to_test_pert/', 'dir')

    mkdir('meshes_to_test_pert')

end

for n = 1:length(name_test)

    wr = ['/home/eit/Mesh_convergence/', name_test{n}, 'img.mat'];

    copyfile(wr, 'meshes_to_test_pert')

end

%% Simulation reconstruction error

near_coord = [0.15, 0.1, 0.13];
far_coord = [0.10,0.17,0.130];
med_coord = [0.13,0.15,0.130];
med2_coord = [0.13,0.14,0.140];
near2_coord = [0.14 0.11 0.12];
near3_coord = [0.14 0.12 0.14]; 
near4_coord = [0.15 0.11 0.14]; 
near5_coord = [0.14 0.10 0.12]; 

pert_test = {near_coord,far_coord,med_coord, med2_coord, near2_coord, near3_coord, near4_coord, near5_coord};
pert_name = {'near', 'far', 'med1', 'med2', 'near2', 'near3', 'near4', 'near5'};

load('names_to_test.mat')

if ~exist('tests_recon', 'dir')

    mkdir('tests_recon')
end


for n = 25% 1:30%length(name_test)

    clear img imgi near_elec Sigma J X c

    cd('/home/yuval/Documents/Mesh_convergence_analysis')

    load([name_test{n},'img.mat'])

    disp(['*************************** Calculating Mesh Number ', num2str(n), ' - ', name_test{n}, ' ******************************'])

    c = find_element_centres(img.fwd_model);

    J = calc_jacobian(img);

    imgh = img;
       
    v1 = fwd_solve(imgh);

    [Mesh_hex,J_hex] = convert_fine2coarse(imgh.fwd_model.elems,imgh.fwd_model.nodes,J,4 * 1e-3); % make element size 0.5mm

    cd('tests_recon')

    %load([name_test{n},'recon.mat'])

    for p = 1:length(pert_test)
        
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
 

        writeVTKcell_orig(['recon_mesh_', num2str(n), '_',pert_name{p}],imgr.fwd_model.elems,imgr.fwd_model.nodes,imgr.elem_data)

        writeVTKcell_orig(['simulated_mesh_', num2str(n), '_',pert_name{p}],imgi.fwd_model.elems,imgi.fwd_model.nodes,imgi.elem_data)
    
    end

    save([name_test{n},'recon.mat'], "pert", '-v7.3')

    close all 

end


%% Calculate WSV 

load('/home/yuval/Documents/Mesh_convergence_analysis/names_to_test.mat')

WSV_true = zeros(length(name_test), length(pert_test)); 
WSV_recon = zeros(length(name_test), length(pert_test)); 

for n = 1:30 %length(name_test) 

    load(['/home/yuval/Documents/Mesh_convergence_analysis/', name_test{n}, 'img.mat'])

    load(['/home/yuval/Documents/Mesh_convergence_analysis/tests_recon/', name_test{n}, 'recon.mat'])

    for p = 1:length(pert_test)

        true_com = pert_test{p}; % true pert centre of mass 
    
        Sigma = pert.Sigma{p}; 
    
        Sigma(isnan(Sigma)) = 0; 
    
        [WSV_true(n,p), WSV_recon(n,p)] = calc_wsv(img, pert.imgi{p},true_com, Sigma, pert.mesh_hex{1}); 

    end 

    disp(['******* Calclated WSV for Mesh Numver ', num2str(n), ' ***********'])

end 

save('WSV', 'WSV_recon', 'WSV_true', 'pert_test', 'pert_name')

%% scale the WSV2 by WSV1 (true pert) 

reg_WSV = WSV_recon./WSV_true; 

min_WSV = reg_WSV(1:26, :) - 1; 

min_WSV = sum(min_WSV, 2); 

plot(reg_WSV) % number closest to 1 is best 
legend(pert_name)


%% Getting the thresholds for image reconstruction visualisation 

for n = 1:30 

    load([name_test{n}, 'recon.mat'])

    for p = 1:length(pert_test)

        pert_threshold(n,p) = max(pert.Sigma{p}); 

    end

    disp(['finished mesh ', num2str(n)])

end 

%% Calc reconstruction 

cd('/home/eit/Mesh_convergence')

for n = 1:length(name_test)

    clear imgr near_elec Sigma J X c

    load([name_test{n},'recon.mat'])

    disp(['*************************** Calculating Mesh ', name_test{n}, ' ******************************'])
    
    imgr = near_elec.imgi;
    imgr.elem_data = near_elec.Sigma;


    figure;
    show_slices(imgr,[inf inf 0.130]);
    title('Regularised Image')
    figure;
    show_slices(near_elec.imgi,[inf inf 0.130]);
    title('Unregularised Image')

end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drafts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forwards per param 

figure
for p = 6:-1:1

    nexttile 

    mat = []; xdat = []; 

    for t = 1:3

        e = params{p}(t); 

        ind = find(ind_map(:,p,:,:,:,:) == e);
        param_combi{p,t}= combi_mesh_fwd(ind,:); 

         % Append the data to mat and xdat
       mat = [mat; param_combi{p, t}(:)];
       xdat = [xdat; t * ones(length(param_combi{p, t}(:)), 1)];
    
       % if isempty(param_combi{p,t})
       % 
       %     continue
       % 
       % else 
       % 
       % 
       %     mat{t} = [param_combi{p,t}(:)];
       %
       % end


    end

    %  violin(mat)

    boxchart(xdat,mat)
    title(param_titles{p})


    xticks(1:3)
    xticklabels(params{p})


end




% %% Forward per param
%
% figure
% 
% for p = 1:length(params)
% 
% 
%     if p == 1
% 
%         param_combi = squeeze(combi_op_fwd(1,:,:,:,:,:,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(2,:,:,:,:,:,:));
% 
%         param_combi2 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(3,:,:,:,:,:,:));
% 
%         param_combi3 = param_combi(:);
% 
%     elseif p == 2
% 
% 
%         param_combi = squeeze(combi_op_fwd(:,1,:,:,:,:,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,2,:,:,:,:,:));
% 
%         param_combi2 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,3,:,:,:,:,:));
% 
%         param_combi3 = param_combi(:);
% 
%     elseif p == 3
% 
%         param_combi = squeeze(combi_op_fwd(:,:,1,:,:,:,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,2,:,:,:,:));
% 
%         param_combi2 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,3,:,:,:,:));
% 
%         param_combi3 = param_combi(:);
% 
%     elseif p == 4
% 
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,1,:,:,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,2,:,:,:));
% 
%         param_combi2 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,3,:,:,:));
% 
%         param_combi3 = param_combi(:);
% 
% 
%     elseif p == 5
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,:,1,:,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,:,2,:,:));
% 
%         param_combi2 = param_combi(:);
% 
%         % param_combi = squeeze(combi_op_fwd(:,:,:,:,3,:,:));
% 
%         param_combi3 = param_combi(:) .* 0;
% 
%     elseif p == 6
% 
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,:,:,1,:));
% 
%         param_combi1 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,:,:,2,:));
% 
%         param_combi2 = param_combi(:);
% 
%         param_combi = squeeze(combi_op_fwd(:,:,:,:,:,3,:));
% 
%         param_combi3 = param_combi(:);
%     end
% 
%     nexttile
% 
%     mm = [param_combi1, param_combi2, param_combi3];
% 
%     boxchart(mm)
% 
%     title(param_titles{p})
% 
%     xticklabels(params{p})
% end




%% Calculate reference forward mesh

name = ['E05_mesh_',num2str(coarse(2), '%.2f'),'_coarse_',num2str(fine(2), '%.2f'),'_fine_',num2str(radius(2)),...
    '_rad_',num2str(facet_ratio(2)),'_ratio_',num2str(facet_angle(2)),'_angle_',num2str(facet_distance(2)),'_dist_'];


cd(inpath)

v1 = []; img = [];

Tetra = readmatrix([name,'tetra.csv']);

Nodes = readmatrix([name,'vertices.csv']);

%Create model

fmdl = mk_common_model('a2C');
fmdl = fmdl.fwd_model;

% add mesh elements and nodes based on MRI

Tetra_orig = Tetra;
Tetra(Tetra(:,5) == 1,:) = [];

% Find boundary and generate eidors model
fmdl.elems = Tetra(:,1:4);
fmdl.nodes = Nodes(:,1:3);
[fmdl.boundary]= dubs3_2(Tetra(:,1:4)); %Tetra(:,1:4)); %find_boundary(Tetra(:,1:4));  %

sig = zeros(length(Tetra),1);
sig(Tetra(:,5) == 2) = 1.79; % csf
sig(Tetra(:,5) == 3) = 0.3; % grey matter
sig(Tetra(:,5) == 4) = 0.15; % white matter


% Add electrodes

cd(outpath)
fmdl.electrode = struct('nodes',[],'z_contact',[]);

img = mk_image(fmdl,sig);

img = add_electrode_yuval(img,elecs,0.00155,0.0002);

[img.fwd_model, unused] = remove_disconnected_nodes(img.fwd_model);
img.fwd_model = remove_unused_nodes(img.fwd_model);

% Add stimulation protocol

img.fwd_model.stimulation = stim;

% Forward Solution

v1 = fwd_solve(img);


% Save forward and Jacobian matrices

ref_fwd = v1.meas;

save([name,'img.mat'],'img');

cd(outpath)
load('All_forwards')

save('All_forwards','Mesh_op_fwd', 'ref_fwd');


%% Comparing forward solutions

load([outpath,'/coarse_forwards_10.mat'])
load([outpath,'/All_forwards_10%.mat'])

coarse_op_fwd_or = coarse_op_fwd;

coarse_op_fwd = cat(1,squeeze(Mesh_op_fwd(1,1,:))', coarse_op_fwd_or);

coarse_op_fwd = coarse_op_fwd([1, 3:end],:);

t_val = [0.45,0.75:0.1:1.85];

for t = 2:size(coarse_op_fwd,1)
    
    c_fine = squeeze(coarse_op_fwd(t-1,:));
    c_coarse = squeeze(coarse_op_fwd(t,:));
    
    E(t-1) = sqrt(sum((c_fine - c_coarse).^2)/1220);
    
end

plot(E, 'marker','o')
xticklabels(t_val)
xticks(1:length(t_val))

%%

cd(outpath)
load('All_forwards.mat')

figure

for p = 1:length(params)
    
    for t = 2
        
        vm = squeeze(abs(Mesh_op_fwd(p,t,:)));
        
        Ev = max((vm - mean(vm))/ mean(vm));
        
        %vr = abs(ref_fwd);
        vr = squeeze(abs(Mesh_op_fwd(p,1,:)));
        
        Ec = max((vm - mean(vr))/ mean(vr));
        
        convg.raw(p,t,:) = [Ev; Ec];
        
        if Ec < Ev
            convg.bin(p,t) = 1;
        else
            convg.bin(p,t) = 0;
        end
        
        
        plot(vm)
        hold on
        plot(vr)
        hold on
    end
end


imagesc(convg.bin)





