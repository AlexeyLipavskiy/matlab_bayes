clc
% close all
clear all
%%
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman');

global days_a_month years sec_in_day mesh_square f_k_flipped s_k mask_square lon_my_mesh...
    lat_my_mesh mesh_b_x_start mesh_b_x_stop mesh_b_y_start mesh_b_y_stop path_to_folder
%% User parameters
years = 2015:2100;


variable_name = 'tos';
path_to_folder = '../CMIP_6/';

exel_list_name = 'list+slp+tos.xls';
save_flag = false;
% save_flag = true;
%% constants
sec_in_day = 60*60*24;
days_a_month = [31,28.25,31,30,31,30,31,31,30,31,30,31];

%% Mask

load rivers_data_year/nor-20_pacif_mask_0.5.mat
%%
% imagesc(lon_full, (lat_full), f_k_north_pacific')
% set(gca,'YDir','normal');
% borders
imagesc(f_k_north_pacific)
%%
path_list = ["../CMIP_6/ssp585/TaiESM1/tos_Omon_TaiESM1_ssp585_r1i1p1f1_gn_208201-208212.nc";
    "../CMIP_6/ssp126/FIO-ESM-2-0/tos_Omon_FIO-ESM-2-0_ssp126_r1i1p1f1_gn_201501-210012.nc";
    "../CMIP_6/ssp126/KACE-1-0-G/tos_Omon_KACE-1-0-G_ssp126_r1i1p1f1_gr_201501-210012.nc";
%     "../Raw/tos_Omon_CESM2-WACCM_ssp126_r1i1p1f1_gr_201501-210012.nc";
    "../CMIP_6/historical/CMCC-CM2-SR5/tos_Omon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412.nc";
    "../CMIP_6/historical/INM-CM5-0/tos_Omon_INM-CM5-0_historical_r1i1p1f1_gr1_185001-189912.nc";
    "../CMIP_6/historical/MPI-ESM1-2-HR/tos_Omon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_185501-185912.nc";
%     "../Raw/tos_Omon_KIOST-ESM_ssp126_r1i1p1f1_gr1_201501-210012.nc"
    
    
];
%%

% path= '..\CMIP_6\ssp585\TaiESM1\tos_Omon_TaiESM1_ssp585_r1i1p1f1_gn_208201-208212.nc';

% path= '..\CMIP_6\ssp126\FIO-ESM-2-0\tos_Omon_FIO-ESM-2-0_ssp126_r1i1p1f1_gn_201501-210012.nc';
path= path_list(2);

a = ncread(path, 'tos');
lat = ncread(path, 'latitude');
lon = ncread(path, 'longitude');
% lat = ncread(path, 'lat');
% lon = ncread(path, 'lon');
tmp = a(:,:,1);
%%
mesh(lat);
figure;
mesh(lon);
%%
path2 = '../CMIP_6/ssp585/TaiESM1/mrro_Lmon_TaiESM1_ssp585_r1i1p1f1_gn_201501-210012.nc';
% path2 = '../CMIP_6/ssp126/MPI-ESM1-2-HR/mrro_Lmon_MPI-ESM1-2-HR_ssp126_r1i1p1f1_gn_203501-203912.nc';

a2 = ncread(path2, 'mrro');
lat2 = ncread(path2, 'lat');
lon2 = ncread(path2, 'lon');
tmp2 = a2(:,:,6);

%%
lon_corr = lon(:,size(lon,1)/2);
lat_corr = lat(round(size(lon,1)/1.01),:);

zer = find(diff(lon_corr) < 0);
if zer < length(lon_corr)/2
    cp = zer + length(lon_corr)/2;
    lon_out = [lon_corr(cp:end-1);lon_corr(1:cp)];
    lon_out(1:length(lon_corr)/2 ) = lon_out(1:length(lon_corr)/2)-360;
    tmp_out = cat(1, tmp(cp:end-1,:,:),tmp(1:cp,:,:));
    
elseif zer > length(lon_corr)/2    
    cp = zer - length(lon_corr)/2;
    lon_out = [lon_corr(cp:end-1);lon_corr(1:cp)]-180; 
    tmp_out = cat(1, tmp(cp:end-1,:,:),tmp(1:cp,:,:));
end








%%
% imagesc(lon_out, lat_corr, tmp_out');
% set(gca,'YDir','normal');
% %%
% 
% imagesc(mesh_b_x, flip(mesh_b_y), f_k_north_pacific);
% set(gca,'YDir','normal');
imagesc(tmp_out);
%%

% [lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(mesh_b_x, mesh_b_y); 

[lon_corr_grid,lat_corr_grid] = ndgrid(lon_out,lat_corr); 
[lon_full_grid,lat_full_grid] = ndgrid(lon_full,lat_full); 


% int_obj = griddedInterpolant(ndgrid(lon_corr,lat_corr), tmp_out);
int_obj = griddedInterpolant(lon_corr_grid, lat_corr_grid, tmp_out);

% f_k_north_pacific_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
tmp_int(:,:) = int_obj(lon_full_grid,lat_full_grid);

%%

imagesc(lon_full, lat_full, (tmp_int.*f_k_north_pacific)');
set(gca,'YDir','normal');
borders













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
mesh_step = 0.5;

x_start_m = min(river_mask(:,1));
x_stop_m = max(river_mask(:,1));
y_start_m = min(river_mask(:,2));
y_stop_m = max(river_mask(:,2));

mesh_b_x_start = fix(x_start_m);
if x_start_m - mesh_b_x_start >= 0.5
    mesh_b_x_start = mesh_b_x_start + 0.5;
end

mesh_b_y_start = fix(y_start_m);
if y_start_m - mesh_b_y_start >= 0.5
    mesh_b_y_start = mesh_b_y_start + 0.5;
end

mesh_b_x_stop = ceil(x_stop_m);
if mesh_b_x_stop - x_stop_m >= 0.5
    mesh_b_x_stop = mesh_b_x_stop - 0.5;
end

mesh_b_y_stop = ceil(y_stop_m);
if mesh_b_y_stop - y_stop_m >= 0.5
    mesh_b_y_stop = mesh_b_y_stop - 0.5;
end


mesh_b_x = mesh_b_x_start : mesh_step : mesh_b_x_stop;
mesh_b_y = mesh_b_y_start : mesh_step : mesh_b_y_stop;

% f_k = zeros(numel(mesh_b_y)-1,numel(mesh_b_x)-1);
parall_degree = pi/180*6377020*cos(deg2rad(mesh_b_y));
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

mesh_lat = mesh_b_x_start + mesh_step/2 : mesh_step : mesh_b_x_stop - mesh_step/2;
mesh_lon = mesh_b_y_start + mesh_step/2 : mesh_step : mesh_b_y_stop - mesh_step/2;

[lon_my_mesh,lat_my_mesh] = ndgrid(mesh_lat,mesh_lon);                                                         % my new mesh 

% clear mesh_mask_intersect mesh_pgon selenga_pgon

%% test part
% figure;
% xa = [96.74 113.25];
% ya = [46.25 53.25]; 
% imagesc(xa,ya,f_k_flipped);
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)

%% CMIP_6

%% read xls list with model names

list_tmp_file = readcell(exel_list_name);
list_tmp_file = list_tmp_file(3:end,:);

% list_hist_models = list_tmp_file(:,1);
% list_hist_marks = list_tmp_file(:,2);
% 
% list_ssp126_models = list_tmp_file(:,13);
% list_ssp126_marks = list_tmp_file(:,14);
% 
% list_ssp245_models = list_tmp_file(:,7);
% list_ssp245_marks = list_tmp_file(:,8);
% 
% list_ssp585_models = list_tmp_file(:,10);
% list_ssp585_marks = list_tmp_file(:,11);

%% preallocation
% var_126 = zeros();


%%
%     SSP_126
list_ssp126_models = list_tmp_file(:,7);
list_ssp126_marks = list_tmp_file(:,8);
[list_of_models_126, var_126] = calc_var(list_ssp126_models, list_ssp126_marks, '/ssp126/', variable_name);
% 
%     SSP_245
list_ssp245_models = list_tmp_file(:,7);
list_ssp245_marks = list_tmp_file(:,8);
[list_of_models_245, var_245] = calc_var(list_ssp245_models, list_ssp245_marks, '/ssp245/', variable_name);

%     SSP_585
list_ssp585_models = list_tmp_file(:,10);
list_ssp585_marks = list_tmp_file(:,11);
[list_of_models_585, var_585] = calc_var(list_ssp585_models, list_ssp585_marks, '/ssp585/', variable_name); 

%%
years_ssp = years;
Rs_126_21 = var_126;
Rs_245_21 = var_245;
Rs_585_21 = var_585;
%%
% save_flag = true;
if save_flag == true
    save("rivers_data_month\"+variable_name+"_"+river_name+"_"+years(1)+"-"+years(end)+"_month_25.01.22.mat",'years_ssp','Rs_126_21'...
    ,'list_of_models_126','Rs_245_21','list_of_models_245','Rs_585_21','list_of_models_585');
end


%%     hist
years = 1979:2014;

list_hist_models = list_tmp_file(:,1);
list_hist_marks = list_tmp_file(:,2);
[list_of_models_hist, var_hist] = calc_var(list_hist_models, list_hist_marks, '/historical/', variable_name); 
%%
years_hist = years;
Rs_hist = var_hist;
%%
% save_flag = true;
if save_flag == true
    save("rivers_data_month\"+variable_name+"_"+river_name+"_"+years(1)+"-"+years(end)+"_month_25.01.22.mat",'years_hist','Rs_hist'...
    ,'list_of_models_hist');
end


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

path = '..\CMIP_6\ssp126\CMCC-CM2-SR5\tos_Omon_CMCC-CM2-SR5_ssp126_r1i1p1f1_gn_201501-210012.nc';
a = ncread(path,'tos');
lon_from_file = ncread(path,'longitude');                                                     % read longitude, maybe not necessary for every file
lat_from_file = ncread(path,'latitude'); 
%%
lon = lon_from_file(:,1);
lat = lat_from_file(1,:);
imagesc(lon, lat,a(:,:,1)')
borders
set(gca,'YDir','normal');


%%
% if save == true
%     save(variable_name+"_"+river_name+"_"+years(1)+"-"+years(end)+"_month_25.01.22.mat",'years_21','P_126_21','R_126_21','Rs_126_21'...
%     ,'list_of_models_126','P_245_21','R_245_21','Rs_245_21','list_of_models_245','P_585_21','R_585_21','Rs_585_21','list_of_models_585');
% end


%%


function[list, var_out] = calc_var(list_models, list_marks, path, var_name)
    global path_to_folder
    count_for_suit_files = 0;
    for count_list = 1:size(list_models,1)
        if string(list_marks(count_list)) == '+'
            count_for_suit_files = count_for_suit_files + 1
%             disp(fullfile(path_to_folder,'/ssp126/',char(list_ssp126_models(count_list))));
            var_out(count_for_suit_files,:) = calculate_cmip6(fullfile(path_to_folder,path,char(list_models(count_list))),var_name);
            list(count_for_suit_files) = list_models(count_list);
        end
    end  
end

function [output] = calculate_cmip6(path,var)
%UNTITLED Summary of this function goes here
% Функция считывает все файлы в папке, на которую указывает путь, фильтрует
% из них те, которые с нужной переменной и подходят по годам (заданным в
% начале скрипта в глобалах). Далее работает с каждым файлом - расчитывает
% средние годовые значения и склеивает их по всем доступным файлам в папке.
% global days_a_month years sec_in_day f_k_flipped s_k mask_square lon_my_mesh lat_my_mesh
global years mask_square


% length_of_var = numel(var);
% list_of_files_tmp = ls(path);
list_of_files_tmp = dir(path);
count_of_files = 0; % counter for suitable files in the folder
disp(path);% to check in command window
% var_year_tmp = 0; % temporary year summ
% year_ind = 1;
month_ind_for_calc_days = 1;
month_total_ind = 1;
for iterator_pre = 3 : size(list_of_files_tmp,1) % start from 3 bc ls give 2 'empty' values in the begining

    file_name_split = strsplit(list_of_files_tmp(iterator_pre).name,'_');% split name for blocks
    file_var = string(file_name_split(1));
    file_years = char(file_name_split(7));
    file_year_start = str2double(file_years(1:4));
    file_year_end = str2double(file_years(8:11));
    
    if file_var == var && file_year_end >= years(1) && file_year_start <= years(end)%find files with matching variable and years

        count_of_files = count_of_files + 1;
%         disp(list_of_files_tmp(iterator_pre,:));
        path_tmp = fullfile(path,list_of_files_tmp(iterator_pre).name); % full path to matched file
        var_tmp = ncread(path_tmp,var);                                                          % read cmip6 file
        lon_from_file = ncread(path_tmp,'lon');                                                     % read longitude, maybe not necessary for every file
        lat_from_file = ncread(path_tmp,'lat'); 

        if count_of_files == 1 % getting years 
            year_start = file_year_start;
        end
        year_stop = file_year_end;
        
        [lon_cmip6,lat_cmip6,lon_ind_cmip6,lat_ind_cmip6] = find_cut_points(lon_from_file,lat_from_file);% different files have diff mesh
        for month_count = 1 : size(var_tmp,3)                                                                         % number of years in file mb different
                                                                                                                  % cycle for every month in file 
            var_month(:,:) = var_tmp(:,:,month_count);                                                                % use one month
            var_month_sum = cut_and_interpolate(var_month,lon_ind_cmip6,lat_ind_cmip6,lon_cmip6,lat_cmip6,month_ind_for_calc_days);
            tmp(month_total_ind) = var_month_sum/mask_square;
            month_total_ind = month_total_ind + 1;

            if mod(month_count,12) == 0                                                                               % when year is full:
                month_ind = 1;
            end  
        end

        if str2double(file_years(5:6)) ~= 1 || str2double(file_years(12:13)) ~= 12% check for errors
            disp('Problem with months. File:   ---------------------------------------------------------------------------------------------------'); 
            disp(list_of_files_tmp(iterator_pre).name);
            count_of_files = 0;
            break;
        end
%         disp('done');
%     elseif file_var ~= var
%         disp('No such variable');
%     elseif file_year_end < years(1)
%         disp('The file is out of years range (hist)');
%     elseif file_year_start > years(end)
%         disp('The file is out of years range (ssp)');
    end
end

if count_of_files == 0
    disp('There are no matching files---------------------------------------------------------------------------------------------------');
%     var_years = NaN;
    output = NaN;
%     error_flag = 1;
else
%     var_years = year_start:year_stop;
    disp(var);
    disp('done');
%     error_flag = 0;
end
var_years = year_start:year_stop;
if years(1) <= 2014
    
    first_ind = find((var_years == years(1))*12);
    output = tmp((first_ind-1)*12 + 1:end);
else
    
end
end

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
