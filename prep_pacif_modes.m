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

% sst_test_ds = deseason(sst,'monthly','dim',3);
% plot(sst_ds)
% hold on
% plot(sst_test_ds)
%%

%%

sst_ds_dt = detrend3(sst_ds,sst_t);
%%














load rivers_data_year\north_pacific_mask.mat
%%
f_k_north_pacific(find(f_k_north_pacific >= 0.5)) = 1;
f_k_north_pacific(find(f_k_north_pacific < 0.5)) = 0;
%%
imagesc((f_k_north_pacific));
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
cor2 = corrcoef(pc_masked(1,:)',ind(:,end));
cor_test = corr(ind, ind);

%%
figure
plot(pc_masked(1,:)/max(pc_masked(1,:)))
hold on
plot(-ind(:,1)/max(ind(:,1)))
%%
% cf = ifft(fft(pc_masked(1,:)).*conj(fft(ind(:,11)))');
% plot(abs(cf))




















