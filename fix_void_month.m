% if ssp data starts from feb, like in NorESM and TaiESM data, the averaging is
% done and calculated data are added to original file

close all
clc
clear all
%%
path_to_hist = '..\CMIP_6\historical\TaiESM1\';
path_to_ssp = '..\CMIP_6\ssp585\TaiESM1\';
var = 'mrros';
averaging_years = 10;

%%

list_of_files_hist = ls(path_to_hist);
list_of_files_ssp = ls(path_to_ssp);
list_of_files_hist = list_of_files_hist(3:end,:);
list_of_files_ssp = list_of_files_ssp(3:end,:);



% paths_ssp = fullfile(path_to_ssp,list_of_files_ssp);

count_pre = 0;
month_count = 0;
left_border = 2014 - averaging_years + 1;

for c = 1 : size(list_of_files_hist,1)
    path_hist = fullfile(path_to_hist,list_of_files_hist(c,:));
    file_name_split = strsplit(list_of_files_hist(c,:),'_');
    file_var = string(file_name_split(1));
    file_years = char(file_name_split(7));
    file_year_start = str2double(file_years(1:4));
    file_year_end = str2double(file_years(8:11));
    num_of_years = file_year_end - file_year_start;
    if file_var == var && file_year_end >= left_border
%         disp(path_hist);
%         count_pre = count_pre + 1;
        
        used_years = file_year_end - max(left_border, file_year_start) + 1;
        
        var_tmp = ncread(path_hist,var);  
        var_hist(:,:,month_count*12 + 1:(month_count + used_years)*12) = var_tmp(:,:,end - (used_years)*12 + 1:end);
        
        month_count = month_count + used_years;
    end
    
end

%%

count_pre = 0;
month_count = 0;
right_border = 2015 + averaging_years;

for c = 1 : size(list_of_files_ssp,1)
    path_ssp = fullfile(path_to_ssp,list_of_files_ssp(c,:));
    file_name_split = strsplit(list_of_files_ssp(c,:),'_');
    file_var = string(file_name_split(1));
    file_years = char(file_name_split(7));
    file_year_start = str2double(file_years(1:4));
    file_year_end = str2double(file_years(8:11));
    num_of_years = file_year_end - file_year_start;
    if file_var == var && file_year_start <= right_border
%         disp(path_ssp);
        count_pre = count_pre + 1;
        
        used_years = min(file_year_end, right_border) - file_year_start + 1;
        
        var_tmp = ncread(path_ssp,var);  
        
        if count_pre == 1 % first file (~corrupted) without January
            num_of_months_corr = min(131, size(var_tmp,3));
            var_ssp(:,:,1:num_of_months_corr) = var_tmp(:,:,1:num_of_months_corr);
            month_count = num_of_months_corr;
            continue
        end
        
        var_ssp(:,:,month_count + 1:month_count + used_years*12) = var_tmp(:,:,1 : (used_years)*12);
    
        month_count = month_count + used_years*12;
    end
    
end

%%
var_full(:,:,1:120) = var_hist;
var_full(:,:,121:240) = var_ssp(:,:,12:end);

inds = 0:averaging_years*2-1;
inds = inds*12 + 1;
var_jans(:,:,1:averaging_years*2) = var_full(:,:,inds);
%%

var_jans_sum = sum(var_jans,3);
var_jans_avg = var_jans_sum./numel(inds);
%%


for c = 1 : size(list_of_files_ssp,1)
    path_ssp = fullfile(path_to_ssp,list_of_files_ssp(c,:));
    file_name_split = strsplit(list_of_files_ssp(c,:),'_');
    file_var = string(file_name_split(1));
    file_years = char(file_name_split(7));
    start_month = file_years(6);
%     file_year_start = str2double(file_years(1:4));
%     file_year_end = str2double(file_years(8:11));
%     num_of_years = file_year_end - file_year_start;
    if file_var == var && start_month == '2'
%         disp(path_ssp);
        var_tmp_corrup = ncread(path_ssp,var);  
        var_name_corrup = list_of_files_ssp(c,:);
    end
    
end
%%
var_fixed(:,:,1) = var_jans_avg;
var_fixed(:,:,2:size(var_tmp_corrup,3)+1) = var_tmp_corrup;

% var_name_fixed = var_name_corrup;
% var_name_fixed(end - 10) = '1';
var_name_fixed = strrep(var_name_corrup,'201502','201501');

%%

path_c = fullfile(path_to_ssp,var_name_corrup);
copyfile(path_c);
ncwrite(var_name_corrup,var,var_fixed);
movefile( var_name_corrup, var_name_fixed)



%%


%%

%%
var_tmp3 = ncread(fullfile('C:\Users\Alex\Desktop\Lab\matlab2020\',var_name_fixed),var);  
% var_tmp2 = ncread(fullfile('C:\Users\Alex\Desktop\Lab\matlab2020\',var_name_corrup),var);  























%%
longitude = ncread(path_ssp,'lon');
latitude = ncread(path_ssp,'lat');

%%
% figure;
% % imagesc(longitude,latitude',(p'),[-50 50]);
% % imagesc(longitude,latitude,var_jans(:,:,1)',[-50 50]);
% % imagesc(longitude,latitude,var_jans(:,:,1)');
% imagesc(longitude,latitude,var_jans(:,:,9)');
% colormap jet;
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)
% 
% figure;
% imagesc(longitude,latitude,var_jans_avg');
% 
% colormap jet;
% set(gca,'YDir','normal'); % change axis direction to normal (default is reverse)

%%

% Plot_cmip6('..\CMIP_6\historical\NorESM2-LM\mrro_Lmon_NorESM2-LM_historical_r1i1p1f1_gn_185001-185912.nc',8);
%%



























