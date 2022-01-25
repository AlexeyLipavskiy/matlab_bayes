clc
% close all
clear all


%%
global days_a_month years sec_in_day mesh_square f_k_flipped s_k mask_square lon_my_mesh lat_my_mesh mesh_b_x_start mesh_b_x_stop mesh_b_y_start mesh_b_y_stop

years = 2016:2100;

%% constants
sec_in_day = 60*60*24;
days_a_month = [31,30,31,30,31,30,31,31,30,31,30,31];
%% Mask
path_mask = '..\r_Selenga.cno';
selenga_mask = dlmread(path_mask);
%% test part
% figure;
% plot(selenga_mask(:,1),selenga_mask(:,2),'.'); % -> 46 : 53.5 lat  
% grid on;                                        %   96.5 : 113.5 long
% figure;
% geoplot(selenga_mask(:,2), selenga_mask(:,1));
% geobasemap colorterrain
%% Mesh (creating regular mesh)
mesh_step = 0.5;
parall_degree = [77465,76761,76057,75342,74627,73900,73173,72435,71697,70948,70199,69439,68679,67909,67138,66358]; % lenght of longitude degree from 46 to 53.5

%bounds
mesh_b_x_start = 96.5;
mesh_b_x_stop = 113.5;
mesh_b_y_start = 46;
mesh_b_y_stop = 53.5;

mesh_b_x = mesh_b_x_start : mesh_step : mesh_b_x_stop;
mesh_b_y = mesh_b_y_start : mesh_step : mesh_b_y_stop;

f_k = zeros(numel(mesh_b_y)-1,numel(mesh_b_x)-1);

%% Counting f_k

selenga_pgon = polyshape(selenga_mask(:,1), selenga_mask(:,2));

for y_count = 1 : numel(mesh_b_y)-1
    
    y1 = (mesh_b_y_stop-(y_count)*mesh_step);
    y2 = y1;
    y3 = (mesh_b_y_stop-(y_count-1)*mesh_step);
    y4 = y3;
        
    for x_count = 1 : numel(mesh_b_x)-1
        x1 = (mesh_b_x_start+(x_count-1)*mesh_step);
        x2 = (mesh_b_x_start+(x_count)*mesh_step);
        x3 = x2;
        x4 = x1;
                 
        mesh_pgon(y_count,x_count) = polyshape([x1 x2 x3 x4],[y1 y2 y3 y4]);
        mesh_square(y_count,x_count) = (mesh_step^2)*(111134.861111*(parall_degree(end-(y_count))+parall_degree(end-(y_count-1)))/2);
                   
    end
end

% figure;
% plot(mesh_pgon);
% hold on;
% plot(selenga_pgon);
% axis image
%% calculating coef. s_k
mesh_mask_intersect = intersect(mesh_pgon,selenga_pgon);
f_k = area(mesh_mask_intersect)./0.25;
f_k_flipped = flip(f_k);
s_k = flip(mesh_square).*ones(numel(mesh_b_y)-1,numel(mesh_b_x)-1).*f_k_flipped;
mask_square = sum(s_k, 'all');
%% mesh

mesh_lat = mesh_b_x_start + mesh_step/2 : mesh_step : mesh_b_x_stop - mesh_step/2;
mesh_lon = mesh_b_y_start + mesh_step/2 : mesh_step : mesh_b_y_stop - mesh_step/2;

[lon_my_mesh,lat_my_mesh] = ndgrid(mesh_lat,mesh_lon);                                                         % my new mesh 

clear mesh_mask_intersect mesh_pgon selenga_pgon

%% test part
% figure;
% xa = [96.74 113.25];
% ya = [46.25 53.25]; 
% imagesc(xa,ya,f_k_flipped);
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)
%% GPCP_2.3

% https://www.ncei.noaa.gov/thredds/catalog/cdr/gpcp_final_agg/catalog.html?dataset=cdr/gpcp_final_agg/CGPC_Final_Aggregation_best.ncd
% https://www1.ncdc.noaa.gov/pub/data/sds/cdr/CDRs/Precipitation_GPCP-Monthly/AlgorithmDescription_01B-34.pdf
% https://www.ncdc.noaa.gov/cdr/atmospheric/precipitation-gpcp-monthly

%%

% path_to_gpcp_folder = '..\GPCP_2.3_data\';
% 
% year_counter = 1;
% gpcp_p_tmp = zeros(144,72);
% gpcp_p_cut = zeros(8,4,12);
% gpcp_p_cut_int = zeros(34,15,12);
% gpcp_p_month_sum = zeros(12,1);
% gpcp_p_year_sum = zeros(size(years,2),1);
% 
% [lon_gpcp,lat_gpcp] = ndgrid(96.25:2.5:113.75,46.25:2.5:53.75);                                                  % part of GPCP2.3 mesh of selenga area
% 
% for year_ind_mrros = years
%     list_of_files_tmp = ls(fullfile(path_to_gpcp_folder,num2str(year_ind_mrros)));
%     list_of_paths_tmp = [path_to_gpcp_folder,num2str(year_ind_mrros),'\'] + string(list_of_files_tmp(3:end,:));        % ls returns two empty fields in the begining
%     
%     for month_ind = 1 : 12     
%         gpcp_p_tmp(:,:) = ncread(list_of_paths_tmp(month_ind,:),'precip');                                       % mm/day       
%         gpcp_p_cut(:,:,month_ind) = gpcp_p_tmp(39:46,55:58);        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %cut required par
%         int_obj = griddedInterpolant(lon_gpcp,lat_gpcp,gpcp_p_cut(:,:,month_ind));                               % create interpolate object   
%         gpcp_p_cut_int(:,:,month_ind) = int_obj(lon_my_mesh,lat_my_mesh);                                       % cutted part interpolated on my new mesh               
%         gpcp_p_month = gpcp_p_cut_int(:,:,month_ind) .* f_k_flipped' .* s_k';                                    % precips in the mask
%         gpcp_p_month_sum(month_ind) = sum(gpcp_p_month,'all')*days_a_month(month_ind);                           % summ of precips for a month            
%     end
%     gpcp_p_year_sum(year_counter) = sum(gpcp_p_month_sum)/mask_square;                                           % summ of all months divided by square
%     
%     year_counter = year_counter + 1;
% end
% 
% disp('GPCP 2.3 done');
%% test part
% P_119(1,:) = gpcp_p_year_sum;
% figure;
% plot(years,gpcp_p_year_sum)

% figure;
% xa = [96.74 113.25];
% ya = [46.25 53.25]; 
% imagesc(xa,ya,gpcp_p_month);
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)
%% netCDF files info

% returns the number of dimensions, variables, and global attributes in a netCDF file
% The function also returns the ID of the dimension defined with unlimited length, if one exists. 
% [ndims, nvars, natts, dimm] = netcdf.inq(ncid_gpc_23);

% varname = netcdf.inqVar(ncid,7); - to recieve the name of a variable
% A_number = netcdf.getVar(ncid,6); - to recieve the variable itself
% ncdisd('path'); - to recieve info about nc file

%% Plot of GPC_2.3 pr. data
% figure;

% % imagesc(longitude,latitude',(p'),[-50 50]);
% imagesc(longitude,latitude,gpcp_p_tmp',[-50 50]);
% colormap jet;
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)

%% EXEL Selenga
% path_xls = '..\rivers';
% selenga_xls = xlsread(path_xls,'selenga');
% count_xls = 1;
% selenga_xls_r = NaN * ones(size(years,2),1);
% for year_ind_xls = years
%     
%     selenga_xls_r(count_xls) = sum(selenga_xls(47 + count_xls,3:14) .* sec_in_day .* (days_a_month));
%     count_xls = count_xls + 1;
% end
% disp('xls data done');
%%
% figure;
% plot(years,selenga_xls_r);


%% CMIP_6

    
% посмотреть данные за декабрь 2020 в гпсп2.3

%% test part
% [hist_pr_years,hist_pr,er] = calculate_cmip6('..\CMIP_6\ssp585\EC-Earth3\','mrros');
% hist_pr = hist_pr/mask_square;
% % [hist_mrro_years,hist_mrro] = calculate_cmip6('..\CMIP_6\historical\AWI-CM-1-1-MR\','mrro');
% % [hist_mrros_years,hist_mrros] = calculate_cmip6('..\CMIP_6\historical\AWI-CM-1-1-MR\','mrros');
% figure;
% plot(hist_pr_years,hist_pr);
% disp(er);
%% read xls list

list_tmp_file = readcell('list.xls');
list_tmp_file = list_tmp_file(3:end,:);

list_hist_models = list_tmp_file(:,1);
list_hist_marks = list_tmp_file(:,2);

list_ssp126_models = list_tmp_file(:,13);
list_ssp126_marks = list_tmp_file(:,14);

list_ssp245_models = list_tmp_file(:,7);
list_ssp245_marks = list_tmp_file(:,8);

list_ssp585_models = list_tmp_file(:,10);
list_ssp585_marks = list_tmp_file(:,11);

% list_of_models_126 = mat2cell('Observed data',1);
% list_of_models_245 = mat2cell('Observed data',1);
% list_of_models_585 = mat2cell('Observed data',1);
%%

%% 
% R_126(1,:) = selenga_xls_r'*1000;
% Rs_126(1,:) = R_126(1,:);
% P_126(1,:) = gpcp_p_year_sum;
% 
% R_245(1,:) = selenga_xls_r'*1000;
% Rs_245(1,:) = R_245(1,:);
% P_245(1,:) = gpcp_p_year_sum;
% 
% R_585(1,:) = selenga_xls_r'*1000;
% Rs_585(1,:) = R_245(1,:);
% P_585(1,:) = gpcp_p_year_sum;

%% SSP_126
intern_count = 0;
for count_list_126 = 1:size(list_tmp_file,1)
    if string(list_ssp126_marks(count_list_126)) == '+'
        intern_count = intern_count + 1
        
        [ssp126_pr_years,ssp126_pr,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp126\',char(list_ssp126_models(count_list_126))),'pr');
        if err_s == 1
            P_126(intern_count,:) = NaN;        
            err_s
        else
            P_126(intern_count,:) = ssp126_pr/mask_square;
        end
     
        [ssp126_mrro_years,ssp126_mrro,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp126\',char(list_ssp126_models(count_list_126))),'mrro');
        if  err_s == 1
            R_126(intern_count,:) = NaN;        
            err_s
        else
            R_126(intern_count,:) = ssp126_mrro/1e12;  
        end
    
        [ssp126_mrros_years,ssp126_mrros,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp126\',char(list_ssp126_models(count_list_126))),'mrros');
        if  err_s == 1
            Rs_126(intern_count,:) = NaN;
            err_s
        else       
            Rs_126(intern_count,:) = ssp126_mrros/1e12;         
        end
        
        list_of_models_126(intern_count) = list_ssp126_models(count_list_126);
    end
end
figure;
plot(years,R_126);

%% SSP_245
intern_count = 0;
for count_list_245 = 1:size(list_tmp_file,1)
    if string(list_ssp245_marks(count_list_245)) == '+'
        intern_count = intern_count + 1
        
        [ssp245_pr_years,ssp245_pr,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp245\',char(list_ssp245_models(count_list_245))),'pr');
        if err_s == 1
            P_245(intern_count,:) = NaN;        
            err_s
        else
            P_245(intern_count,:) = ssp245_pr/mask_square;
        end
     
        [ssp245_mrro_years,ssp245_mrro,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp245\',char(list_ssp245_models(count_list_245))),'mrro');
        if  err_s == 1
            R_245(intern_count,:) = NaN;        
            err_s
        else
            R_245(intern_count,:) = ssp245_mrro/1e12;  
        end
    
        [ssp245_mrros_years,ssp245_mrros,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp245\',char(list_ssp245_models(count_list_245))),'mrros');
        if  err_s == 1
            Rs_245(intern_count,:) = NaN;
            err_s
        else       
            Rs_245(intern_count,:) = ssp245_mrros/1e12;         
        end
        
        list_of_models_245(intern_count) = list_ssp245_models(count_list_245);
    end
end
figure;
plot(years,R_245);



%% SSP_585

intern_count = 0;
for count_list_585 = 1:size(list_tmp_file,1)
    if string(list_ssp585_marks(count_list_585)) == '+'
        intern_count = intern_count + 1
        
        [ssp585_pr_years,ssp585_pr,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp585\',char(list_ssp585_models(count_list_585))),'pr');
        if err_s == 1
            P_585(intern_count,:) = NaN;        
            err_s
        else
            P_585(intern_count,:) = ssp585_pr/mask_square;
        end
     
        [ssp585_mrro_years,ssp585_mrro,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp585\',char(list_ssp585_models(count_list_585))),'mrro');
        if  err_s == 1
            R_585(intern_count,:) = NaN;        
            err_s
        else
            R_585(intern_count,:) = ssp585_mrro/1e12;  
        end
    
        [ssp585_mrros_years,ssp585_mrros,err_s] = calculate_cmip6(fullfile('..\CMIP_6\ssp585\',char(list_ssp585_models(count_list_585))),'mrros');
        if  err_s == 1
            Rs_585(intern_count,:) = NaN;
            err_s
        else       
            Rs_585(intern_count,:) = ssp585_mrros/1e12;         
        end
        
        list_of_models_585(intern_count) = list_ssp585_models(count_list_585);
    end
end
figure;
plot(years,R_585);


%%



years_21 = years;
P_126_21 = P_126;
R_126_21 = R_126;
Rs_126_21 = Rs_126;
P_245_21 = P_245;
R_245_21 = R_245;
Rs_245_21 = Rs_245;
P_585_21 = P_585;
R_585_21 = R_585;
Rs_585_21 = Rs_585;

% save('ssp126+245+585_V2(2016-2100).mat','years_21','P_126_21','R_126_21','Rs_126_21','list_of_models_126','P_245_21','R_245_21','Rs_245_21'...
%     ,'list_of_models_245','P_585_21','R_585_21','Rs_585_21','list_of_models_585');



%%
function [var_years,output,error_flag] = calculate_cmip6(path,var)
%UNTITLED Summary of this function goes here
% Функция считывает все файлы в папке, на которую указывает путь, фильтрует
% из них те, которые с нужной переменной и подходят по годам (заданным в
% начале скрипта в глобалах). Далее работает с каждым файлом - расчитывает
% средние годовые значения и склеивает их по всем доступным файлам в папке.
% global days_a_month years sec_in_day f_k_flipped s_k mask_square lon_my_mesh lat_my_mesh
global years

% length_of_var = numel(var);
list_of_files_tmp = ls(path);
count_pre = 0; % counter for suitable files in the folder
disp(path);% to check in command window
var_year_tmp = 0; % temporary year summ
year_ind = 1;
month_ind = 1;

for iterator_pre = 3 : size(list_of_files_tmp,1) % start from 3 bc ls give 2 empty values in the begining
    file_name_split = strsplit(list_of_files_tmp(iterator_pre,:),'_');% split name for blocks 
    file_var = string(file_name_split(1));
    file_years = char(file_name_split(7));
    file_year_start = str2double(file_years(1:4));
    file_year_end = str2double(file_years(8:11));
    if file_var == var && file_year_end >= years(1) && file_year_start <= years(end)%find files with matching variable and years
        count_pre = count_pre + 1;
%         disp(list_of_files_tmp(iterator_pre,:));
        path_tmp = fullfile(path,list_of_files_tmp(iterator_pre,:)); % full path to matched file
        var_tmp = ncread(path_tmp,var);                                                          % read cmip6 file
        lon_from_file = ncread(path_tmp,'lon');                                                     % read longitude, maybe not necessary for every file
        lat_from_file = ncread(path_tmp,'lat'); 

        if count_pre == 1 % getting years 
            year_start = years(1);
        end
        year_stop = years(end);
        
        [lon_cmip6,lat_cmip6,lon_ind_cmip6,lat_ind_cmip6] = find_cut_points(lon_from_file,lat_from_file);% different files have diff mesh
        if file_year_start <= years(1)
            month_start = (years(1) - file_year_start)*12;
        else
            month_start =1;
        end
         if file_year_end >= years(end)
             month_stop = (years(end) - file_year_start +1 )*12;
         else
             month_stop = size(var_tmp,3);
         end
        
        
        for month_count = month_start +1 : month_stop                                                                         % number of years in file mb different
%                      month_count                                                                % cycle for every month in file 
            var_month(:,:) = var_tmp(:,:,month_count);                                                                % use one month
            var_month_sum = cut_and_interpolate(var_month,lon_ind_cmip6,lat_ind_cmip6,lon_cmip6,lat_cmip6,month_ind);
            var_year_tmp = var_year_tmp + var_month_sum;                                                              % summ of flow for a year
            if mod(month_count,12) == 0                                                                               % when year is full:
                output(year_ind) = var_year_tmp;                                                                      % get result
                var_year_tmp = 0;
                year_ind = year_ind +1;
                month_ind = 1;
            end  
        end

        if str2double(file_years(5:6)) ~= 1 || str2double(file_years(12:13)) ~= 12% check for errors
            disp('Problem with months. File:   ---------------------------------------------------------------------------------------------------'); 
            disp(list_of_files_tmp(iterator_pre,:));
            count_pre = 0;
            break;
        end
%         disp('done');
%     elseif file_var ~= var
%         disp('Other variable');
%     elseif file_year_end < years(1)
%         disp('The file is out of years range (hist)');
%     elseif file_year_start > years(end)
%         disp('The file is out of years range (ssp)');
    end
end

if count_pre == 0
    disp('There are no matching files');
    var_years = NaN;
    output = NaN;
    error_flag = 1;
else
    var_years = year_start:year_stop;
    disp(var);
    disp('done');
    error_flag = 0;
end
% figure;
% plot(var_years,output);
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
int_obj = griddedInterpolant(lon_cmip6,lat_cmip6,var_month_cut(:,:));                             % create interpolate object   
var_month_cut_int(:,:) = int_obj(lon_my_mesh,lat_my_mesh);                                        % cutted part interpolated on my new mesh
var_month_cut_int = var_month_cut_int .* f_k_flipped' .* s_k';                                  % flow in the mask
var_month_sum = sum(var_month_cut_int,'all')*days_a_month(month_ind)*sec_in_day;                % summ of flow for a month
end


% function [fixed] = fix_years(first_years,first,second_years,second)
% % Cut required part of data according to years used in analys
% 
% global years
% %days_a_month sec_in_day f_k_flipped s_k lon_my_mesh lat_my_mesh
%   
% first_ind = find(first_years == years(1));
% first_part = first(first_ind:end);
% 
% second_ind_start = find(second_years == first_years(end) + 1);
% second_ind_end = find(second_years == years(end));
% second_part = second(second_ind_start:second_ind_end);
% 
% fixed = [first_part,second_part];
% end
