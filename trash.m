% imagesc(lon_full, (lat_full), f_k_north_pacific')
% set(gca,'YDir','normal');
% borders
% imagesc(f_k_north_pacific)
%%
% path_list = ["../CMIP_6/historical/BCC-CSM2-MR/tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412.nc";
%     "../Raw/tas_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_185001-201412.nc";
%     "../Raw/tas_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412.nc";
%     "../Raw/tas_Amon_FIO-ESM-2-0_historical_r1i1p1f1_gn_185001-201412.nc";
%     "../Raw/tas_Amon_INM-CM5-0_historical_r1i1p1f1_gr1_195001-201412.nc";
%     "../Raw/tas_Amon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185001-185412.nc";
%     "../Raw/tas_Amon_TaiESM1_historical_r1i1p1f1_gn_185001-201412.nc";
% %     "../Raw/tos_Omon_KIOST-ESM_ssp126_r1i1p1f1_gr1_201501-210012.nc"
%     
%     
% ];
%%

% % path= '..\CMIP_6\ssp585\TaiESM1\tos_Omon_TaiESM1_ssp585_r1i1p1f1_gn_208201-208212.nc';
% 
% path= '..\CMIP_6\ssp126\FIO-ESM-2-0\tos_Omon_FIO-ESM-2-0_ssp126_r1i1p1f1_gn_201501-210012.nc';
% path= path_list(1);
% % 
% a = ncread(path, 'tas');
% % % lat = ncread(path, 'latitude');
% % % lon = ncread(path, 'longitude');
% lat = ncread(path, 'lat');
% lon = ncread(path, 'lon');
% tmp = a(:,:,1);
%%
% mesh(lat);
% figure;
% mesh(lon);
%%
% path2 = '../CMIP_6/ssp585/TaiESM1/mrro_Lmon_TaiESM1_ssp585_r1i1p1f1_gn_201501-210012.nc';
% % path2 = '../CMIP_6/ssp126/MPI-ESM1-2-HR/mrro_Lmon_MPI-ESM1-2-HR_ssp126_r1i1p1f1_gn_203501-203912.nc';
% 
% a2 = ncread(path2, 'mrro');
% lat2 = ncread(path2, 'lat');
% lon2 = ncread(path2, 'lon');
% tmp2 = a2(:,:,6);
% %%
% imagesc(lon, lat, tmp');
% borders
% set(gca,'YDir','normal');
% figure;
% imagesc(lon2, lat2, tmp2);
% borders


%%
% lat_mask = lat_full;
% lon_mask = linspace(0, 359.5, 720);
% %%
% save rivers_data_year/nor-20_pacif_mask_0.5_shift.mat f_k_north_pacific s_k_north_pacific lon_mask lat_mask
% %%
% s_k_north_pacific = cat(1, s_k_north_pacific(361:end,:), s_k_north_pacific(1:360,:));
% %%
% f_k_north_pacific = cat(1, f_k_north_pacific(361:end,:), f_k_north_pacific(1:360,:));
%%
% imagesc(lon_mask, lat_full, (f_k_north_pacific)');
% set(gca,'YDir','normal');
% borders
% figure;
% imagesc(lon, lat, tmp');
% set(gca,'YDir','normal');
% borders
% % lon_corr = lon(:,size(lon,1)/2);
% % % lat_corr = lat(round(size(lon,1)/1.01),:);
% % lat_corr = flip(lat(201,:));
%% zero meridian shift

% zer = find(diff(lon_corr) < 0);
% if zer < length(lon_corr)/2
%     cp = zer + length(lon_corr)/2;
%     lon_out = [lon_corr(cp+1:end);lon_corr(1:cp)];
%     lon_out(1:length(lon_corr)/2 ) = lon_out(1:length(lon_corr)/2)-360;
%     tmp_out = cat(1, tmp(cp:end-1,:,:),tmp(1:cp,:,:));
%     
% elseif zer > length(lon_corr)/2    
%     cp = zer - length(lon_corr)/2;
%     lon_out = [lon_corr(cp+1:end);lon_corr(1:cp)]; 
%     lon_out(1:ceil(length(lon_corr)/2)) = lon_out(1:ceil(length(lon_corr)/2))-360;
%     tmp_out = cat(1, tmp(cp:end-1,:,:),tmp(1:cp,:,:));
% end


%%
% imagesc(lon_out, lat_corr, tmp_out');
% set(gca,'YDir','normal');
% %%
% 
% imagesc(mesh_b_x, flip(mesh_b_y), f_k_north_pacific);
% set(gca,'YDir','normal');
% imagesc(tmp_out);

%%
% [lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(mesh_b_x, mesh_b_y); 

% [lon_mask_grid,lat_mask_grid] = ndgrid(lon_mask,lat_mask); 
% [lon_var_grid,lat_var_grid] = ndgrid(lon,lat); 
% 
% 
% % int_obj = griddedInterpolant(ndgrid(lon_corr,lat_corr), tmp_out);
% int_obj = griddedInterpolant(lon_mask_grid, lat_mask_grid, f_k_north_pacific);
% 
% % f_k_north_pacific_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
% mask_int(:,:) = int_obj(lon_var_grid,lat_var_grid);
% 
% %%
% figure
% % imagesc(lon_full, lat_full, (tmp_int.*f_k_north_pacific)');
% 
% imagesc(lon, lat, (logical(mask_int).*tmp)');
% 
% set(gca,'YDir','normal');
% borders


%%
% tmp2(tmp2 == 0) = NaN;
% 
% % imagesc(b, lat(1,:), tmp')
% % set(gca,'YDir','normal');
% % figure;
% % imagesc(lon2, lat2, tmp2')
% % set(gca,'YDir','normal');
% imagesc( tmp_out');
% set(gca,'YDir','normal');
% figure;
% 
% imagesc(tmp2' );
% set(gca,'YDir','normal');
% 
% %%
% b = lon(:,size(lon,1)/2);
% 
% %%
% mesh(b, lat(1,:), tmp');
% figure;
% mesh(lon);



%% Mesh (creating regular mesh)

%bounds
% mesh_step = 0.5;
% 
% x_start_m = min(river_mask(:,1));
% x_stop_m = max(river_mask(:,1));
% y_start_m = min(river_mask(:,2));
% y_stop_m = max(river_mask(:,2));
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
% % f_k = zeros(numel(mesh_b_y)-1,numel(mesh_b_x)-1);
% parall_degree = pi/180*6377020*cos(deg2rad(mesh_b_y));
%% Counting f_k

% river_pgon = polyshape(river_mask(:,1), river_mask(:,2));
% 
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
% end
%%
% figure;
% plot(mesh_pgon);
% hold on;
% plot(selenga_pgon);
% axis image
% grid on; 
% xlabel('Latitude');
% ylabel('Longitude');
%% calculating coef. s_k
% mesh_mask_intersect = intersect(mesh_pgon,river_pgon);
% f_k = area(mesh_mask_intersect)./0.25;
% f_k_flipped = flip(f_k);
% s_k = flip(mesh_square).*ones(numel(mesh_b_y)-1,numel(mesh_b_x)-1).*f_k_flipped;
% mask_square = sum(s_k, 'all');
%% mesh

% mesh_lat = mesh_b_x_start + mesh_step/2 : mesh_step : mesh_b_x_stop - mesh_step/2;
% mesh_lon = mesh_b_y_start + mesh_step/2 : mesh_step : mesh_b_y_stop - mesh_step/2;
% 
% [lon_my_mesh,lat_my_mesh] = ndgrid(mesh_lat,mesh_lon);                                                         % my new mesh 

% clear mesh_mask_intersect mesh_pgon selenga_pgon

%% test part
% figure;
% xa = [96.74 113.25];
% ya = [46.25 53.25]; 
% imagesc(xa,ya,f_k_flipped);
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)



%%
% year_ind = 1;
% var_year_tmp = 0;
% month_count = 1;
% for i = 1: size(var_126,2)
%     
%             var_year_tmp = var_year_tmp + var_126(6, i);                                                              % summ of flow for a year
%             if mod(month_count,12) == 0                                                                               % when year is full:
%                 output2(year_ind) = var_year_tmp;                                                                      % get result
%                 var_year_tmp = 0;
%                 year_ind = year_ind +1;
% %                 month_ind = 1;
%             end 
%             month_count = month_count +1;
% end
% z = 2015:2100;
% plot(z, output2)
%%
% clear all;
% %%
a = ncinfo('../Merra_2_data/MERRA2_100.instM_2d_asm_Nx.198001.nc4');
% 
% b = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','SLP');
% lat = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lat');
% lon = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lon');
% %%
% imagesc(lat, lon, b')
% 
% %%
% % mask_square = 1;
% merra2_years = 1980:2014;
% path_to_merra2_folder = '../Merra2/';
% 
% year_counter = 0;
% merra2_psl_tmp = zeros(576, 361);
% merra2_psl_month_sum = zeros(12*size(merra2_years,2),1);
% lat_from_file = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lat');
% lon_from_file = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lon');
% ls_tmp = dir(path_to_merra2_folder);
% 
% %%
% for year_ind = merra2_years
%     
%     
%     for month_ind = 1 : 12     
%         path_merra2_tmp = strcat(path_to_merra2_folder, ls_tmp(year_counter*12 + month_ind + 2).name);
%         merra2_psl_tmp(:,:) = ncread(path_merra2_tmp,'SLP');                                                         % mm/day  
% 
%         [lon_gpcp,lat_gpcp,lon_ind_gpcp,lat_ind_gpcp] = find_cut_points(...
%             lon_from_file,lat_from_file);  % different files have diff mesh
%  
%         merra2_psl_month_sum(12*year_counter + month_ind) = ...
%             cut_and_interpolate(merra2_psl_tmp,...
%             lon_ind_gpcp,lat_ind_gpcp,lon_gpcp,lat_gpcp,month_ind)/...
%             sec_in_day/mask_square;
%         
%     end
%     year_counter = year_counter + 1;
% end
% % overall units are mm/year
% disp('pr observed data (GPCP 2.3) done');
%%
% clc
% path = '/home/alex/Downloads/Merra_1/MERRA2_100.instM_2d_lfo_Nx.198001.nc4';
% path = '/home/alex/Downloads/Merra_2/MERRA2_100.instU_3d_ana_Np.198001.nc4';
% path = '/home/alex/Downloads/Merra_3/MERRA2_100.instU_2d_lfo_Nx.198001.nc4';
% path = '/home/alex/Downloads/Merra_4/MERRA2_100.instM_3d_asm_Np.198001.nc4';
% path = '/home/alex/Downloads/Merra_5/MERRA2_100.instU_3d_asm_Np.198001.nc4';
% path = '/home/alex/Downloads/Merra_6/MERRA2_100.tavgM_2d_slv_Nx.198001.nc4';
% path = '/home/alex/Downloads/Merra_7/MERRA2_100.tavgU_2d_slv_Nx.198001.nc4';
path = '../Merra_2_data/MERRA2_100.instM_3d_ana_Np.198001.nc4';
% path = '/home/alex/Downloads/Merra_9/MERRA2.tavgC_3d_ltm_Np.198101_201001.nc4';



ncdisp(path);
%%
aa = ncread(path, 'lev');
bb = ncread(path, 'lon');
cc = ncread(path, 'lat');














%%
clear all
close all
clc
%%

load index.mat
%%
% ind_t_tmp = sst_t(1309:1740);
% R_des = deseason(R_hist_m(1,:), ind_t_tmp);

ind_cut = ind(349:780,:)';

%%
R_y = month_to_year_data(R_hist_m);
ind_y = month_to_year_data(ind_cut);



%%
[cor,p_val] = corr(R_hist_m(1:12,:)',ind_cut');
[cor_1,p_val_1] = corr(R_y',ind_y');
cor2 = corrcoef(R_hist_m(1,:),ind_cut(3,:));
cor_test = corr(ind_cut, ind_cut);

%%

[1 "North Atlantic Oscillation (NAO)";
 2   "East Atlantic (EA)";
 3   "West Pacific (WP)";
 4   "East Pacific-North Pacific (EP-NP)";
 5   "Pacific/North American (PNA)";
 6   "East Atlantic/Western Russia";
 7   "Scandinavia (SCAND)";
 8   "Tropical/Northern Hemisphere (TNH)";
 9   "Polar/Eurasia";
 10   "Pacific Transition (PT)";
 11   "Pacific Decadal Oscillation (PDO)";
 12   "Oceanic Niño Index (ONI)"]












function [lon_cut,lat_cut,lon_cut_ind,lat_cut_ind] = find_cut_points(given_lon,given_lat)
%   Input vars are lon and lat of a given file. Using coordinates of
%   borders of my_grid (defined by mask) cut_points are found to then cut from
%   given file. Cut points have to be outside the my_mesh area to perform
%   interpolation
global mesh_b_x_start mesh_b_x_stop mesh_b_y_start mesh_b_y_stop
lon_start = mesh_b_x_start;
lon_stop = mesh_b_x_stop;
lat_start = mesh_b_y_start;
lat_stop = mesh_b_y_stop;

given_mesh_step_lon =  given_lon(2) - given_lon(1);
given_mesh_step_lat =  given_lat(2) - given_lat(1);

lon_cut_start_ind = find(given_lon < lon_start,1,'last');
lon_cut_stop_ind = find(given_lon > lon_stop,1,'first');
lat_cut_start_ind = find(given_lat < lat_start,1,'last');
lat_cut_stop_ind = find(given_lat > lat_stop,1,'first');

lon_cut_start_value = given_lon(lon_cut_start_ind);
lon_cut_stop_value = given_lon(lon_cut_stop_ind);
lat_cut_start_value = given_lat(lat_cut_start_ind);
lat_cut_stop_value = given_lat(lat_cut_stop_ind);

[lon_cut,lat_cut] = ndgrid(lon_cut_start_value:given_mesh_step_lon:lon_cut_stop_value,lat_cut_start_value:given_mesh_step_lat:lat_cut_stop_value+0.001);

lon_cut_ind = lon_cut_start_ind : lon_cut_stop_ind;
lat_cut_ind = lat_cut_start_ind : lat_cut_stop_ind;
end


function [var_month_sum] = cut_and_interpolate(var_month,lon_ind_cmip6,lat_ind_cmip6,lon_cmip6,lat_cmip6,month_ind)
global days_a_month sec_in_day f_k_flipped s_k lon_my_mesh lat_my_mesh
  
var_month_cut(:,:) = var_month(lon_ind_cmip6,lat_ind_cmip6);                                    % cut required part

nans = isnan(var_month_cut);                                                                    % find nans in file and replace it with a 0
var_month_cut(nans == 1) = 0;

int_obj = griddedInterpolant(lon_cmip6,lat_cmip6,var_month_cut(:,:));                             % create interpolate object   
var_month_cut_int(:,:) = int_obj(lon_my_mesh,lat_my_mesh);                                        % cutted part interpolated on my new mesh
var_month_cut_int = var_month_cut_int .* f_k_flipped' .* s_k';                                  % flow in the mask
var_month_sum = sum(var_month_cut_int,'all')*days_a_month(month_ind)*sec_in_day;                % summ of flow for a month
end


function[output] = month_to_year_data(mon_data)
    if size(mon_data,1) > size(mon_data,2)
        mon_data = mon_data';
    end
    n_of_models = size(mon_data, 1);
    output = zeros(n_of_models, ceil(size(mon_data, 2)/12));
    year_data_tmp = zeros(n_of_models, 1);
    year_ind = 1;
    month_count = 1;
    
    for i = 1: size(mon_data,2)

                year_data_tmp = year_data_tmp + mon_data(:, i);                 % summ of flow for a year
                if mod(month_count,12) == 0                                % when year is full:
                    output(:, year_ind) = year_data_tmp;                     % get result
                    year_data_tmp = zeros(n_of_models, 1);
                    year_ind = year_ind +1;
    %                 month_ind = 1;
                end 
                month_count = month_count +1;
    end
end
