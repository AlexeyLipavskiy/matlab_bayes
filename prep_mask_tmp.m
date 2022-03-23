clc
close all
clear all
%%



load rivers_data_year\north_pacific_mask.mat
%%
f_k_north_pacific(find(f_k_north_pacific >= 0.5)) = 1;
f_k_north_pacific(find(f_k_north_pacific < 0.5)) = 0;
%%
% imagesc((f_k_north_pacific));
%%
f_k_north_pacific = flip(f_k_north_pacific);
%%
na_lat_tmp = 0.5:0.5:67;
na_lon_tmp = -179.5:0.5:180;

% imagesc(na_lon_tmp, na_lat_tmp, (f_k_north_pacific));
% borders
% set(gca,'YDir','normal');
%%
% imagesc(sst_lon,sst_lat,sst_ds_dt(:,:,1)')

% 
% lat_new = 67.5:-1:0.5;
% lon_new = -179.5:1:179.5;
% %%
% [lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 
% % [na_lon_tmp,na_lat_tmp] = ndgrid(mesh_b_x, mesh_b_y); 
% 
% int_obj = griddedInterpolant(na_lon_tmp, na_lat_tmp, f_k_north_pacific');                             % create interpolate object   
% f_k_north_pacific_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
%%








% imagesc(na_lon_tmp, na_lat_tmp, flip(f_k_north_atlantic));
% borders
% set(gca,'YDir','normal');
%%
% figure;
% imagesc(lon_new,lat_new, (f_k_north_pacific_int)');
% borders
% set(gca,'YDir','normal');
%%

% imagesc(lat_new,lon_new, flip(f_k_north_atlantic_int'));
% borders
% set(gca,'YDir','normal');
mask1 = zeros(720, 360);
lon_full = -180:0.5:179.5;
lat_full = -90:0.5:89.5;
%%
% mask1(find(sst_lon == lon_new),find(sst_lat == lat_new)) = f_k_north_atlantic_int;
for n = 1:length(mesh_b_x)
    for t = 1:length(mesh_b_y)
        mask1(lon_full == mesh_b_x(n),lat_full == mesh_b_y(t)) = s_k_north_pacific(t,n);
        
    end

end
%%
imagesc(mask1);
%% 
%зануляем маску от экватора до 20 гр с.ш.
mask1(:,180:220) = 0;

%%
imagesc(lon_full,lat_full, mask1');
borders
set(gca,'YDir','normal');
%%
s_k_north_pacific = mask1;
%%
imagesc(lon_full,lat_full, s_k_north_pacific');
borders
set(gca,'YDir','normal');
figure;
imagesc(lon_full,lat_full, f_k_north_pacific');
borders
set(gca,'YDir','normal');

% mask_na = logical(mask1);
%%
% save rivers_data_year/nor-20_pacif_mask_0.5.mat f_k_north_pacific s_k_north_pacific lat_full lon_full
























