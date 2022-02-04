% clc
% clear all
%%
% path = '..\CMIP_6\ssp585\ACCESS-CM2';
path = '../Raw/';
var = 'sst';
%%
global years  f_k_north_pacific
load rivers_data_year/nor-20_pacif_mask_0.5_shift.mat
% years = 2015:2100;
years = 1979:2014;
%%

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
        var_tmp = ncread(path_tmp,var);                                                          % read cmip6 file
      

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
    output_2(:,:,mon) = output(:,:,mon).*logical(mask_int);
%     fillmissing( output(:,:,mon),'constant',0);
end
%%
% imagesc(lon_from_file,lat_from_file, output(:,:,671)')
% borders
% set(gca,'YDir','normal');
%%
datenum_t = double(time + datenum(years(1),1,0));

output_ds = deseason(output_2,datenum_t);

output_ds_dt = detrend3(output_ds,datenum_t);

[eof_maps_masked,pc_masked,expv_masked] = eof(output_ds_dt,1);
%%

plot(pc_masked)
%%
pdo_tmp = pdo(349:780);
%%
pc_tmp = pc_masked(1:432);
%%
cor2 = corrcoef(pc_tmp',pdo_tmp);
%%
load index.mat
%%
for n =1:size(var_hist,1)
   cor(:,:,n) = corrcoef(pdo_tmp', var_hist(n,:)); 
end

%%
% year_ind = 1;
% var_year_tmp = 0;
% month_count = 1;
% for i = 1: size(output,2)
%     
%             var_year_tmp = var_year_tmp + output(i);                                                              % summ of flow for a year
%             if mod(month_count,12) == 0                                                                               % when year is full:
%                 output3(year_ind) = var_year_tmp;                                                                      % get result
%                 var_year_tmp = 0;
%                 year_ind = year_ind +1;
% %                 month_ind = 1;
%             end 
%             month_count = month_count +1;
% end
%%
% z = 1979:2014;
% plot(years, output3)
%%


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







