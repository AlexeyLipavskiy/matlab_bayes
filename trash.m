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






























