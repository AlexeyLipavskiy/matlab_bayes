% sort .nc files from CMIP6 project from one folder (Raw)
% to folders like
% CMIP_6 -> ssp126 -> GFDL


% clc
close all
clear all

%%
years = 1979:2020;
path_raw = '..\Raw\';
path_dest = '..\CMIP_6\';


%%

list_of_files_raw = ls(path_raw);
list_of_files_raw = list_of_files_raw(3:end,:);
% list_of_folders = ls(path_dest);
% list_of_folders = list_of_folders(3:end,:);


%%

n_of_files = length(list_of_files_raw);
% n_of_files = 10;
f = waitbar(0,'Please wait...');



for count = 1:n_of_files
    name_split = strsplit(list_of_files_raw(count,:),'_');
    exp_name = string(name_split(4));
    path_exp = path_dest+exp_name;
    
    model_name = string(name_split(3));
    list_of_models = string(ls(path_exp));
    
    if isempty(find(list_of_models == model_name,1)) % esli net takoy papki s imenem modeli
        cd (path_exp)
        mkdir (model_name)
        if isfile(fullfile(path_exp,model_name,list_of_files_raw(count,:))) % check if file already exists
            disp(list_of_files_raw(count,:));
            disp('This file already exists in this folder');
            continue
        end
        movefile (fullfile("..\",path_raw,list_of_files_raw(count,:)) , fullfile("..\",path_exp,model_name))
        cd ..\ % move to CMIP_6
    else
        if isfile(fullfile(path_exp,model_name,list_of_files_raw(count,:))) % check if file already exists
            disp(list_of_files_raw(count,:));
            disp('This file already exists in this folder');
            continue
        end

        movefile (fullfile(path_raw,list_of_files_raw(count,:)) , fullfile(path_exp,model_name))       
    end
    waitbar(count/n_of_files,f);
end

close(f)

%%

% ext = '*.nc'; % extension you care about
% d = dir(fullfile(path_raw,ext));
% filenames = {d(:).name}';
% 
% mask1 = ~cellfun('isempty',regexp(filenames,'hist'));
% mask2 = ~cellfun('isempty',regexp(filenames,'type 4'));
% mask3 = ~mask1 & ~mask2;


% src = fullfile(pathname,filenames(mask1));
% dest = fullfile('X',filenames(mask1));
% cellfun(@copyfile,src,dest);
% 
% src = fullfile(pathname,filenames(mask2));
% dest = fullfile('Y',filenames(mask2));
% cellfun(@copyfile,src,dest);
% 
% src = fullfile(pathname,filenames(mask3));
% dest = fullfile('Z',filenames(mask3));
% cellfun(@copyfile,src,dest);









