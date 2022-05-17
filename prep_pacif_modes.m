clc
close all
clear all
%%
load index.mat

%% const
years = 1979:2020;
path_index = '../Tele_index/HadISST_sst.nc';
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

% sst_test_ds = deseason(sst,'monthly','dim',3);
plot(sst_ds)
% hold on
% plot(sst_test_ds)
%%

%%

sst_ds_dt = detrend3(sst_ds,sst_t);
%%














load rivers_data_year/north_atlantic_mask.mat
%%
% f_k_north_pacific(find(f_k_north_atlantic >= 0.5)) = 1;
% f_k_north_pacific(find(f_k_north_atlantic < 0.5)) = 0;
%%
imagesc((f_k_north_atlantic));
%%
f_k_north_pacific = flip(f_k_north_pacific);
%%
na_lat_tmp = 0.5:0.5:67;
na_lon_tmp = -179.5:0.5:180;

% imagesc(na_lon_tmp, na_lat_tmp, (f_k_north_pacific));
% borders
% set(gca,'YDir','normal');
%%
% imagesc(sst_lon,sst_lat,sst_ds_dt(:,:,1)')


lat_new = 67.5:-1:0.5;
lon_new = -179.5:1:179.5;
%%
[lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
[na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(mesh_b_x, mesh_b_y); 

int_obj = griddedInterpolant(na_lon_tmp, na_lat_tmp, f_k_north_pacific');                             % create interpolate object   
f_k_north_pacific_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
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
mask1 = zeros(length(sst_lon), length(sst_lat));
%%
% mask1(find(sst_lon == lon_new),find(sst_lat == lat_new)) = f_k_north_atlantic_int;
for n = 1:length(lon_new)
    for m = 1:length(lat_new)
        mask1(sst_lon == lon_new(n),sst_lat == lat_new(m)) = f_k_north_pacific_int(n,m);
        
    end

end
%%
%% 
%зануляем маску от экватора до 20 гр с.ш.
mask1(:,71:100) = 0;

%%
imagesc(sst_lon,sst_lat, mask1');
borders
set(gca,'YDir','normal');
%%


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
psl_ds_masked = zeros(size(sst_ds));
for v = 1:size(sst_ds,3)
    psl_ds_masked(:,:,v) = sst_ds_dt(:,:,v).*mask1;
end
%%
[eof_maps_masked,pc_masked_tmp,expv_masked] = eof(psl_ds_masked,1);
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
[cor,p_val] = corr(pc_masked_tmp',ind);
cor2 = corrcoef(pc_masked_tmp(1,:)',ind(:,end));
cor_test = corr(ind, ind);

%%
figure
plot(pc_masked_tmp(1,:)/max(pc_masked_tmp(1,:)))
hold on
plot(-ind(:,1)/max(ind(:,1)))
%%
% cf = ifft(fft(pc_masked(1,:)).*conj(fft(ind(:,11)))');
% plot(abs(cf))
%%

plot(sst_dt, pc_masked_tmp)
%%
pc_pdo_cut = pc_masked_tmp(349:780);


%% processing slp data (Merra2)


merra2_years = 1980:2014;
path_to_merra2_folder = '../Merra2/';


psl = zeros(576, 361, 420);
merra2_psl_month_sum = zeros(12*size(merra2_years,2),1);
lat_from_file = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lat');
lon_from_file = ncread('../Merra2/MERRA2_100.instM_2d_asm_Nx.198001.nc4','lon');
ls_tmp = dir(path_to_merra2_folder);
psl_dt = zeros(420,1);
%%
year_counter = 0;
for year_ind = merra2_years
    
    
    for month_ind = 1 : 12     
        path_merra2_tmp = strcat(path_to_merra2_folder, ls_tmp(year_counter*12 + month_ind + 2).name);
        psl(:,:, year_counter*12 + month_ind) = ncread(path_merra2_tmp,'SLP');    
        psl_dt(year_counter*12 + month_ind) = ncread(path_merra2_tmp,'time'); 

        
    end
    year_counter = year_counter + 1;
end
% overall units are mm/year
disp('psl data done');
%%
imagesc(psl(:,:,1));
%%
psl_t = double(psl_dt + datenum(1980,1,0));
psl_dt_norm = datetime(1980,1,0,'Format','yyyy-MM') + caldays(ceil(psl_dt));
%%
psl_ds = deseason(psl, psl_t);
psl_ds_dt = detrend(psl_ds, psl_t);
%%
load rivers_data_year/north_atlantic_mask.mat

% f_k_north_pacific(find(f_k_north_pacific >= 0.5)) = 1;
% f_k_north_pacific(find(f_k_north_pacific < 0.5)) = 0;

% imagesc((f_k_north_pacific));

f_k_north_pacific = flip(f_k_north_atlantic);
%% mesh of the mask
na_lat_tmp = 0:0.5:68.5;
na_lon_tmp = -97.5:0.5:12.5;

% imagesc(na_lon_tmp, na_lat_tmp, (f_k_north_pacific));
% borders
% set(gca,'YDir','normal');
%% mesh of the data (psl here)
% imagesc(sst_lon,sst_lat,sst_ds_dt(:,:,1)')


lat_new = 0:0.5:68.5;
lon_new = -97.5:0.625:12.5;
%%
[lon_new_tmp,lat_new_tmp] = ndgrid(lon_new,lat_new); 
[na_lon_tmp,na_lat_tmp] = ndgrid(na_lon_tmp,na_lat_tmp); 
% [na_lon_tmp,na_lat_tmp] = ndgrid(mesh_b_x, mesh_b_y); 

int_obj = griddedInterpolant(na_lon_tmp, na_lat_tmp, f_k_north_pacific');                             % create interpolate object   
f_k_north_pacific_int(:,:) = int_obj(lon_new_tmp,lat_new_tmp); 
%%
% imagesc(na_lon_tmp, na_lat_tmp, flip(f_k_north_atlantic));
% borders
% set(gca,'YDir','normal');
%%
figure;
imagesc(lon_new,lat_new, (f_k_north_pacific_int)');
borders
set(gca,'YDir','normal');

%%
% mask1 = zeros(length(sst_lon), length(sst_lat));
mask1 = zeros(length(lon_from_file), length(lat_from_file));
%%
% mask1(find(sst_lon == lon_new),find(sst_lat == lat_new)) = f_k_north_atlantic_int;
for n = 1:length(lon_new)
    for m = 1:length(lat_new)
        mask1(lon_from_file == lon_new(n),lat_from_file == lat_new(m)) = f_k_north_pacific_int(n,m);
        
    end

end
%%
%% 
%зануляем маску от экватора до 20 гр с.ш.
mask1(:,170:221) = 0;

%%
imagesc(lon_from_file,lat_from_file, mask1');
borders
set(gca,'YDir','normal');
%%
imagesc(lon_from_file, lat_from_file, psl(:,:,2)')
borders
set(gca,'YDir','normal');

figure;

imagesc(lon_from_file, lat_from_file, sst(:,:,2)')
borders
set(gca,'YDir','normal');
%%
mask_na = logical(mask1);

%%
psl_ds_masked = zeros(size(psl));
for v = 1:size(psl,3)
    psl_ds_masked(:,:,v) = psl(:,:,v).*mask1;
end
%%
imagesc(lon_from_file, lat_from_file, psl_ds_masked(:,:,2)')
borders
set(gca,'YDir','normal');

%%
[eof_maps_masked,pc_masked_tmp,expv_masked] = eof(psl_ds_masked,12);

%%

ind_cut = ind(361:780,:);

%%
[cor,p_val] = corr(pc_masked_tmp',ind_cut);
cor2 = corrcoef(pc_masked_tmp(1,:)',ind_cut(:,2));
cor_test = corr(ind_cut, ind_cut);






























