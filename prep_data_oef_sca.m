clc
% close all
clear all
%%
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman');

% global days_a_month years sec_in_day mesh_square f_k_flipped s_k mask_square lon_my_mesh...
%     lat_my_mesh mesh_b_x_start mesh_b_x_stop mesh_b_y_start mesh_b_y_stop path_to_folder
global path_to_folder years mask_int_obj
%% User parameters
years = 2015:2100;


variable_name = 'zg';
path_to_folder = '../CMIP_6/';

exel_list_name = 'list+slp+tas.xls';
save_flag = false;
% save_flag = true;
%% constants
% sec_in_day = 60*60*24;
% days_a_month = [31,28.25,31,30,31,30,31,31,30,31,30,31];
% sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_187001-202112.nc
%% Mask
% load rivers_data_year/nor-20_atlan_mask_0.5_shift.mat
% load rivers_data_year/nor-20_pacif_mask_0.5_shift.mat
load rivers_data_year/SCA_mask.mat
%%
% imagesc(lon_mask,lat_mask, sca_mask');
% borders
% set(gca,'YDir','normal');

% imagesc(lon_mask,lat_mask, f_k_north_pacific');
% borders
% set(gca,'YDir','normal');
%%
[lon_mask_grid,lat_mask_grid] = ndgrid(lon_mask,lat_mask); 
% % mask_int_obj = griddedInterpolant(lon_mask_grid, lat_mask_grid, f_k_north_pacific);
% mask_int_obj = griddedInterpolant(lon_mask_grid, lat_mask_grid, f_k_north_atlantic);
mask_int_obj = griddedInterpolant(lon_mask_grid, lat_mask_grid, f_k_sca); 

%% CMIP_6

%% read xls list with model names

list_tmp_file = readcell(exel_list_name);
list_tmp_file = list_tmp_file(3:end,:);

list_hist_models = list_tmp_file(:,1);
list_hist_marks = list_tmp_file(:,2);

list_ssp126_models = list_tmp_file(:,13);
list_ssp126_marks = list_tmp_file(:,14);

list_ssp245_models = list_tmp_file(:,7);
list_ssp245_marks = list_tmp_file(:,8);

list_ssp585_models = list_tmp_file(:,10);
list_ssp585_marks = list_tmp_file(:,11);

%% inds to cut
% if we need DJF-mean (for example) 
%using only 1,2,12 months
inds = 3:11;
for n = 1:numel(years)
    ind_to_cut((n-1)*numel(inds) + 1 : n*numel(inds)) = inds + (n-1)*12;
    
end


%%
years = 2015:2100;
% %     SSP_126
list_ssp126_models = list_tmp_file(:,7);
list_ssp126_marks = list_tmp_file(:,8);
[list_of_models_126, var_126] = calc_var(list_ssp126_models, list_ssp126_marks, '/ssp126/', variable_name);

% %     SSP_245
list_ssp245_models = list_tmp_file(:,7);
list_ssp245_marks = list_tmp_file(:,8);
[list_of_models_245, var_245] = calc_var(list_ssp245_models, list_ssp245_marks, '/ssp245/', variable_name);

%%     SSP_585
list_ssp585_models = list_tmp_file(:,10);
list_ssp585_marks = list_tmp_file(:,11);
[list_of_models_585, var_585] = calc_var(list_ssp585_models, list_ssp585_marks, '/ssp585/', variable_name); 

%%
years_ssp = years;
sca_126_21 = var_126;
sca_245_21 = var_245;
sca_585_21 = var_585;
%% modification of nao


sca_126_21_cut = sca_126_21;
sca_245_21_cut = sca_245_21;
sca_585_21_cut = sca_585_21;


sca_126_21_cut(:,ind_to_cut) = 0;
sca_245_21_cut(:,ind_to_cut) = 0;
sca_585_21_cut(:,ind_to_cut) = 0;


%%
% save_flag = true;
if save_flag == true
    save("rivers_data_month/"+"sca_cut("+variable_name+")_"+years(1)+"-"+years(end)+"_month_20.09.22.mat",'years_ssp','sca_126_21_cut'...
    ,'list_of_models_126','sca_245_21_cut','list_of_models_245','sca_585_21_cut','list_of_models_585');
end

if save_flag == true
    save("rivers_data_month/"+"sca("+variable_name+")_"+years(1)+"-"+years(end)+"_month_20.09.22.mat",'years_ssp','sca_126_21'...
    ,'list_of_models_126','sca_245_21','list_of_models_245','sca_585_21','list_of_models_585');
end

disp("------SSP DONE---------------------------------------------------------");
%%     hist
% years = 1979:2014;
% 
% list_hist_models = list_tmp_file(:,1);
% list_hist_marks = list_tmp_file(:,2);
% [list_of_models_hist, var_hist] = calc_var(list_hist_models, list_hist_marks, '/historical/', variable_name); 
% %%
% years_hist = years;
% % pdo_hist = pc_pdo_cut;    КОСТЫЛЬ
% % pdo_hist(2:19,:) = var_hist;
% 
% 
% 
% load index.mat
% nao_hist = sca(349:780)';
% nao_hist(2:19,:) = var_hist;
% 
% %% modification of nao
% 
% nao_hist_cut = nao_hist;
% nao_hist_cut(:,ind_to_cut) = 0;
% plot(nao_hist_cut(1,:))
% hold on;
% plot(nao_hist(1,:))
% %%
% 
% sca_hist = nao_hist;
% sca_hist_cut = nao_hist_cut;
% 
% 
% %%
% % save_flag = true;
% if save_flag == true
%     save("rivers_data_month/"+"sca("+variable_name+")_"+years(1)+"-"+years(end)+"_month_20.09.22.mat",'years_hist','sca_hist'...
%     ,'list_of_models_hist');
% end
% disp("-----------------------DONE --------------------------------------");
% %%
% 
% if save_flag == true
%     save("rivers_data_month/"+"sca_cut("+variable_name+")_"+years(1)+"-"+years(end)+"_month_20.09.22.mat",'years_hist','sca_hist_cut'...
%     ,'list_of_models_hist');
% end
% 















%%

function[list, var_out] = calc_var(list_models, list_marks, path, var_name)
    global path_to_folder
    count_for_suit_files = 0;
    for count_list = 1:size(list_models,1)
        if string(list_marks(count_list)) == '+'
            count_for_suit_files = count_for_suit_files + 1
%             disp(fullfile(path_to_folder,'/ssp126/',char(list_ssp126_models(count_list))));
            var_out(count_for_suit_files,:) = calculate_cmip6_eof(fullfile(path_to_folder,path,char(list_models(count_list))),var_name);
            list(count_for_suit_files) = list_models(count_list);
        end
    end  
end

function [pc_masked] = calculate_cmip6_eof(path,var)
%UNTITLED Summary of this function goes here
% Функция считывает все файлы в папке, на которую указывает путь, фильтрует
% из них те, которые с нужной переменной и подходят по годам (заданным в
% начале скрипта в глобалах). Далее работает с каждым файлом - расчитывает
% средние годовые значения и склеивает их по всем доступным файлам в папке.
% global days_a_month years sec_in_day f_k_flipped s_k mask_square lon_my_mesh lat_my_mesh
global years mask_int_obj


%%

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
        var_tmp_4d = ncread(path_tmp,var);                                                          % read cmip6 file
        var_tmp_4d_shift = permute(var_tmp_4d, [1 2 4 3]);
        
        
        var_tmp = var_tmp_4d_shift(:,:,:,8);

        if count_of_files == 1 % getting years 
            year_start = file_year_start;
            output_tmp = var_tmp;
            lon_from_file = ncread(path_tmp,'lon');                                                     % read longitude, maybe not necessary for every file
            lat_from_file = ncread(path_tmp,'lat');
            time_from_file = ncread(path_tmp,'time');
        else
            output_tmp = cat(3, output_tmp, var_tmp);
            time_from_file = [time_from_file; ncread(path_tmp,'time')];
        end
        year_stop = file_year_end;
        
        
%         [lon_cmip6,lat_cmip6,lon_ind_cmip6,lat_ind_cmip6] = find_cut_points(lon_from_file,lat_from_file);% different files have diff mesh
%         for month_count = 1 : size(var_tmp,3)                                                                         % number of years in file mb different
%                                                                                                                   % cycle for every month in file 
%             var_month(:,:) = var_tmp(:,:,month_count);                                                                % use one month
%             var_month_sum = cut_and_interpolate(var_month,lon_ind_cmip6,lat_ind_cmip6,lon_cmip6,lat_cmip6,month_ind_for_calc_days);
%             tmp(month_total_ind) = var_month_sum/mask_square;
%             month_total_ind = month_total_ind + 1;
% 
%             if mod(month_count,12) == 0                                                                               % when year is full:
%                 month_ind = 1;
%             end  
%         end

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
%%
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
    output = output_tmp(:,:,(first_ind-1)*12 + 1:end);
    time = time_from_file((first_ind-1)*12 + 1:end);
else
    last_ind = find((var_years == years(end)));
    output = output_tmp(:,:,1:last_ind*12);
    time = time_from_file(1:last_ind*12);
end
% figure;
% plot(var_years,output);
%%
 
[lon_var_grid,lat_var_grid] = ndgrid(lon_from_file,lat_from_file); 
mask_int(:,:) = mask_int_obj(lon_var_grid,lat_var_grid);
%%
for mon = 1:size(output,3)
    output(:,:,mon) = output(:,:,mon).*logical(mask_int);
end
%%
% imagesc(lon_from_file,lat_from_file, output(:,:,671)')
% borders
% set(gca,'YDir','normal');
%%
datenum_t = double(time + datenum(years(1),1,0));

output_ds = deseason(output,datenum_t);

output_ds_dt = detrend3(output_ds,datenum_t);

[eof_maps_masked,pc_masked,expv_masked] = eof(output_ds_dt,1);
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