clc
close all
clear all
%%



% load rivers_data_year/north_pacific_mask.mat
load rivers_data_year/north_atlantic_mask.mat
%%
% f_k_north_pacific(find(f_k_north_pacific >= 0.5)) = 1;
% f_k_north_pacific(find(f_k_north_pacific < 0.5)) = 0;
%%
% imagesc((f_k_north_pacific));
%%
% f_k_north_pacific = flip(f_k_north_pacific);
%%
% f_k_north_pacific = mask_na;
s_k_north_atlantic = flip(s_k_north_atlantic);
%%
na_lat_tmp = 68.5:-0.5:0;
na_lon_tmp = -98:0.5:12;

% % na_lat_tmp = 0.5:0.5:67;
% na_lon_tmp = -179.5:0.5:180;
% % na_lon_tmp = 180:-0.5:179.5;
% na_lat_tmp = 67:-0.5:0.5;
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
mask2 = zeros(720, 360);
lon_full = -180:0.5:179.5;
lat_full = -90:0.5:89.5;
%%

% for n = 1:length(mesh_b_x)
%     for t = 1:length(mesh_b_y)
%         mask1(lon_full == mesh_b_x(n),lat_full == mesh_b_y(t)) = f_k_north_pacific(t,n);
%         
%     end
% 
% end

for n = 1:length(na_lon_tmp)
    for t = 1:length(na_lat_tmp)
        mask1(lon_full == na_lon_tmp(n),lat_full == na_lat_tmp(t)) = ...
            f_k_north_atlantic(t,n);
        
    end

end

for n = 1:length(na_lon_tmp)
    for t = 1:length(na_lat_tmp)
        mask2(lon_full == na_lon_tmp(n),lat_full == na_lat_tmp(t)) = ...
            flip(s_k_north_atlantic(t,n));
        
    end

end


%%
imagesc(mask1);
%% 
%зануляем маску от экватора до 20 гр с.ш.
mask1(:,180:220) = 0;
mask2(:,180:220) = 0;
%%
imagesc(lon_full,lat_full, mask1');
borders
set(gca,'YDir','normal');
%%
f_k_north_atlantic = mask1;
s_k_north_atlantic = mask2;
%%
imagesc(lon_full,lat_full, s_k_north_atlantic');
borders
set(gca,'YDir','normal');
figure;
imagesc(lon_full,lat_full, f_k_north_atlantic');
borders
set(gca,'YDir','normal');

%% shifting
tmp = f_k_north_atlantic(1:360,:);
f_k_north_atlantic(1:360,:) = f_k_north_atlantic(361:end,:);
f_k_north_atlantic(361:end,:) = tmp;
%%
tmp = s_k_north_atlantic(1:360,:);
s_k_north_atlantic(1:360,:) = s_k_north_atlantic(361:end,:);
s_k_north_atlantic(361:end,:) = tmp;

%%
lon_full = 0:0.5:359.5;

%%
imagesc(lon_full,lat_full, f_k_north_atlantic');
borders
set(gca,'YDir','normal');
figure;
imagesc(lon_full,lat_full, s_k_north_atlantic');
borders
set(gca,'YDir','normal');
% mask_na = logical(mask1);
%%

lat_mask = lat_full;
lon_mask = lon_full;
%%
save rivers_data_year/nor-20_atlan_mask_0.5_shift.mat ...
    f_k_north_atlantic s_k_north_atlantic lon_mask lat_mask

%% SCA for Merra2


path = '../Merra_2_data/MERRA2_100.instM_3d_ana_Np.198001.nc4';
h = ncread(path, 'H');
gh = h(:,:,1);
lon_merra = ncread(path, 'lon');
lat_merra = ncread(path, 'lat');

%%
imagesc(lon_merra,lat_merra, gh');
borders
set(gca,'YDir','normal');
%%

mask = zeros(576, 361);
mask(358:529, 221:356) = ones(172, 136);
%%
imagesc(lon_merra,lat_merra, (gh.*mask)');
borders
set(gca,'YDir','normal');
%%
sca_mask = mask;
%%
save rivers_data_year/SCA_mask.mat ...
    sca_mask



























































