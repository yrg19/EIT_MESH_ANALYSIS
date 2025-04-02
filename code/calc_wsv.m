%% Calculate Weighted Spatial Variance 

function [WSV_true, WSV_recon] = calc_wsv(orig_img, true_pert_img, true_com, recon_sigma, Mesh_hex)

% outputs are the WSV values for the true perturbation and the
% reconstructed
% inputs are:
% original img no perturbation
% true perturbation image with perturbation
% reconstructed sigma values
% center of mass of true pertrubation coordinate
% Mesh hex from the reconstruction

%% Caclulate WSV for true perturbation

c_tet = find_element_centres(orig_img.fwd_model);

temp1 = true_pert_img.elem_data - orig_img.elem_data; % true pert

V = get_elem_volume(orig_img.fwd_model);

summy = sum(V.*temp1.^2);
w = V.*temp1.^2/summy;
WSV_true = sqrt(sum(w.*((c_tet(:,1)-true_com(1)).^2 + (c_tet(:,2)-true_com(2)).^2 + (c_tet(:,3)-true_com(3)).^2)));


%% Calculate WSV for reconstruction

c_hex = find_element_centres_hex(Mesh_hex.Hex, Mesh_hex.Nodes); % get this for hex

temp2 = recon_sigma;  % reconstruction

%  temp1(abs(img1.elem_data)<max(abs(img1.elem_data))*0.5) = 0;

temp2(abs(recon_sigma)<max(abs(recon_sigma))*0.5) = 0; % play around with absolute vs signed thresholding

V = ones(size(Mesh_hex.Hex)).* (Mesh_hex.d) .^3;


summy = sum(V.*temp2.^2);
w = V.*temp2.^2/summy;
WSV_recon = sqrt(sum(w.*((c_hex(:,1)-true_com(1)).^2 + (c_hex(:,2)-true_com(2)).^2 + (c_hex(:,3)-true_com(3)).^2)));

end

