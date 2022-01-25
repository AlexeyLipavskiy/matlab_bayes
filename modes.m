clc
close all
clear all
%%
load index.mat

%% const
years = 1979:2020;
path_index = '..\Tele_index\HadISST_sst.nc';
global sst_lat sst_lon
%% read data

sst = ncread(path_index, 'sst');
%%
sst_time = ncread(path_index, 'time');
sst_lat = ncread(path_index, 'latitude');
sst_lon = ncread(path_index, 'longitude');


%%
sst_t = double(ncread(path_index,'time') + datenum(1870,1,0));
sst_dt = datetime(1870,1,0,'Format','yyyy-MM') + caldays(ceil(sst_time));

%%
sst_tmp = sst(:,:,1);
sst_tmp(find (sst_tmp < -500)) = NaN;

month_number=1;
%% 
sst(find (sst < -500)) = NaN;
sst(:,:,1:960) =[];

%%
sst_dt(1:960) =[];
sst_t(1:960) =[];
% [eof_maps,pc,expvar] = eof(sst,20);
%%
% plot_sst(eof_maps,1)
% plot_sst(sst_tmp,1)
%%
% imagescn(sst_lon,sst_lat,eof_maps(:,:,1)');
% cmocean('delta','pivot',0)
%%
% plot(sst_dt, pc(1,:))
%%
sst_ds = deseason(sst,sst_t);
%%

sst_ds_dt = detrend3(sst_ds,sst_t);
%%














load north_atlantic_mask.mat
%%
f_k_north_atlantic(find(f_k_north_atlantic >= 0.5)) = 1;
f_k_north_atlantic(find(f_k_north_atlantic < 0.5)) = 0;
% imagesc((f_k_north_atlantic));
%%
na_lat_tmp = 0:0.5:68.5;
na_lon_tmp = -97.5:0.5:12.5;
f_k_north_atlantic = flip(f_k_north_atlantic);
imagesc(na_lon_tmp, na_lat_tmp, (f_k_north_atlantic));
borders
set(gca,'YDir','normal');
%%
% imagesc(sst_lon,sst_lat,sst_ds_dt(:,:,1)')

% lat_new = 0.5:1:68.5;
lat_new = 68.5:-1:0.5;
lon_new = -97.5:1:12.5;
%%
[lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
[na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 

int_obj = griddedInterpolant(na_lon_tmp, na_lat_tmp, f_k_north_atlantic');                             % create interpolate object   
f_k_north_atlantic_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
%%








% imagesc(na_lon_tmp, na_lat_tmp, flip(f_k_north_atlantic));
% borders
% set(gca,'YDir','normal');
%%
figure;
imagesc(lon_new,lat_new, (f_k_north_atlantic_int)');
borders
set(gca,'YDir','normal');
%%

% imagesc(lat_new,lon_new, flip(f_k_north_atlantic_int'));
% borders
% set(gca,'YDir','normal');
mask1 = zeros(length(sst_lon), length(sst_lat));
%%
% mask1(find(sst_lon == lon_new),find(sst_lat == lat_new)) = f_k_north_atlantic_int;
for n = 1:length(lon_new)
    for m = 1:length(lat_new)
        mask1(sst_lon == lon_new(n),sst_lat == lat_new(m)) = f_k_north_atlantic_int(n,m);
        
    end

end
%%

imagesc(sst_lon,sst_lat, mask1');
borders
set(gca,'YDir','normal');
%%
% for n = 1:length(sst_lon)
%     for m = 1:length(sst_lat)
%         if mask1(n,m) == 1 && isnan(sst_tmp(n,m))
%             mask1(n,m) = 0;
%         end
%     end
% end
% for n = 1:length(sst_lon)
%     for m = 1:length(sst_lat)
%         if mask1(n,m) == 0
%             mask1(n,m) = NaN;
%         end
%     end
% end
%%
% imagesc(sst_lon,sst_lat, mask1');
% set(gca,'YDir','normal');
% figure;
% 
% imagesc(sst_lon,sst_lat, sst_tmp');
% set(gca,'YDir','normal');
%%
% % imagesc(sst_lon,sst_lat, mask1.*sst_tmp);
% imagesc(sst_lon,sst_lat,(mask1.*sst_tmp)');
% borders
% set(gca,'YDir','normal');
%%
mask_na = logical(mask1);
%%
% save('north_atlantic_mask.mat','s_k_north_atlantic','f_k_north_atlantic','mask_na');

%%
sst_ds_masked = zeros(size(sst_ds));
for v = 1:size(sst_ds,3)
    sst_ds_masked(:,:,v) = sst_ds_dt(:,:,v).*mask1;
end
%%
[eof_maps_masked,pc_masked,expv_masked] = eof(sst_ds_masked,11);
% [eof_maps,pc,expv] = eof(sst_ds,11);
%%

% Plot the first mode:
% figure
% imagesc(sst_lon,sst_lat,eof_maps(:,:,1)')
% axis xy image
% cmocean('curl','pivot')
% title 'The first EOF mode'
% 
% figure
% imagesc(sst_lon,sst_lat,eof_maps_masked(:,:,1)')
% axis xy image
% cmocean('curl','pivot')
% title 'The first EOF mode'
%%

%%
% nao(end-2:end) = [];
ind(end-2:end,:) = [];
ind(:,12) = [];
%%
[cor,p_val] = corr(pc_masked',ind);
cor2 = corrcoef(pc(1,:)',ind(:,end));
cor_test = corr(ind, ind);
%%
imagesc(cor_test);
%%

cor2 = corrcoef(pc(1,:),pc(1,:));

plot(pc_masked(1,:)/max(pc_masked(1,:)))
hold on
plot(ind(:,1)/max(ind(:,1)))

%%
% % Determine which grid points correspond to land:
% [lat,lon] = cdtgrid([1/2 1/2]);
% land = island(lat,lon);
% 
% % Set land values to NaN:
% z = zeros(size(lat));
% z(land) = 1;
% 
% 
% %%
% 
% z(180:end,:) = 0;
% 
% %%
% 
% % Plot the masked dataset:
% figure
% imagescn(lon,lat,z)
% hold on
% borders
% 
% %%
% ax = worldmap('Pacific');
% % geoshow('landareas.shp')
% % setm(gca,'ffacecolor',rgb('ocean blue'))
% 
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow(ax, land, 'FaceColor', [0.5 0.7 0.8])
% %%
% 
% load coastlines
%%
% path_ocean_borders = '../ocean_borders/goas_v01.shp';
% ocean_shape = shaperead(path_ocean_borders);
% %%
% ocean_shape = ocean_shape(4);
% ocean_mask(:,1) = ocean_shape.X;
% ocean_mask(:,2) = ocean_shape.Y;
% %%
% figure;
% geoplot(ocean_mask(:,2),ocean_mask(:,1));
% % geoplot(north_atlantic_mask(:,2),north_atlantic_mask(:,1),'.');
% geobasemap colorterrain
% % geoshow(north_atlantic_mask(:,2),north_atlantic_mask(:,1));
% 
% %%
% %% Mesh (creating regular mesh)
% 
% %bounds
% mesh_step = 0.5;
% 
% x_start_m = min(ocean_mask(:,1));
% x_stop_m = max(ocean_mask(:,1));
% y_start_m = min(ocean_mask(:,2));
% y_stop_m = max(ocean_mask(:,2));
% 
% mesh_b_x_start = fix(x_start_m);
% if x_start_m - mesh_b_x_start >= 0.5
%     mesh_b_x_start = mesh_b_x_start + 0.5;
% end
% 
% mesh_b_y_start = fix(y_start_m);
% if y_start_m - mesh_b_y_start >= 0.5
%     mesh_b_y_start = mesh_b_y_start + 0.5;
% end
% 
% mesh_b_x_stop = ceil(x_stop_m);
% if mesh_b_x_stop - x_stop_m >= 0.5
%     mesh_b_x_stop = mesh_b_x_stop - 0.5;
% end
% 
% mesh_b_y_stop = ceil(y_stop_m);
% if mesh_b_y_stop - y_stop_m >= 0.5
%     mesh_b_y_stop = mesh_b_y_stop - 0.5;
% end
% 
% 
% mesh_b_x = mesh_b_x_start : mesh_step : mesh_b_x_stop;
% mesh_b_y = mesh_b_y_start : mesh_step : mesh_b_y_stop;
% 
% f_k = zeros(numel(mesh_b_y)-1,numel(mesh_b_x)-1);
% parall_degree = pi/180*6377020*cos(deg2rad(mesh_b_y));
% %% Counting f_k
% 
% river_pgon = polyshape(ocean_mask(:,1), ocean_mask(:,2));
% disp('river_pgon');
% f = waitbar(0);
% for y_count = 1 : numel(mesh_b_y)-1
%     
%     y1 = (mesh_b_y_stop-(y_count)*mesh_step);
%     y2 = y1;
%     y3 = (mesh_b_y_stop-(y_count-1)*mesh_step);
%     y4 = y3;
%         
%     for x_count = 1 : numel(mesh_b_x)-1
%         x1 = (mesh_b_x_start+(x_count-1)*mesh_step);
%         x2 = (mesh_b_x_start+(x_count)*mesh_step);
%         x3 = x2;
%         x4 = x1;
%                  
%         mesh_pgon(y_count,x_count) = polyshape([x1 x2 x3 x4],[y1 y2 y3 y4]);
%         mesh_square(y_count,x_count) = (mesh_step^2)*(111134.861111*(parall_degree(end-(y_count))+parall_degree(end-(y_count-1)))/2);
%                    
%     end
%     waitbar(y_count/numel(mesh_b_y)-1, f);
% end
% close(f);
% disp('f_k');
% % figure;
% % plot(mesh_pgon);
% % hold on;
% % plot(selenga_pgon);
% % axis image
% % grid on; 
% % xlabel('Latitude');
% % ylabel('Longitude');
% 
% mesh_mask_intersect = intersect(mesh_pgon,river_pgon);
% disp('mesh_mask_intersect');
% f_k = area(mesh_mask_intersect)./0.25;
% disp('area');
% f_k_flipped = flip(f_k);
% s_k = flip(mesh_square).*ones(numel(mesh_b_y)-1,numel(mesh_b_x)-1).*f_k_flipped;
% mask_square = sum(s_k, 'all');
% %%
% s_k_north_pacific = s_k;
% f_k_north_pacific = f_k;
% %%
% save('north_pacific_mask.mat','s_k_north_pacific','f_k_north_pacific','mesh_b_x','mesh_b_y');
































function plot_sst(sst,month_number)
% 

global sst_lat sst_lon
imagesc(sst_lon,sst_lat,sst(:,:,month_number)');
colormap jet;
set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)
xlabel('Longitude, \circ');
ylabel('Latitude, \circ');
colorbar
end




