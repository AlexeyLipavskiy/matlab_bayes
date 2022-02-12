% write nc-file in order to test "test function" with HadISST data as it
% was model data
%%
years = 1979:2020;
path_index = '..\Tele_index\HadISST_sst.nc';
global sst_lat sst_lon
%% read data

sst = ncread(path_index, 'sst');
%%
sst_t = double(ncread(path_index,'time') + datenum(1870,1,0));
sst_dt = datetime(1870,1,0,'Format','yyyy-MM') + caldays(ceil(sst_time));
%%

sst_cut = flip(sst(:,:,1309:1740),2);
%%
sst_lat = ncread(path_index, 'latitude');
sst_lon = ncread(path_index, 'longitude');
%%
imagesc(sst_cut(:,:,1));

% set(gca,'YDir','normal');
figure;
imagesc(f_k_north_pacific);
%%
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


for m = 1:size(sst_cut,3)
    ttt = sst_cut(1:size(sst_cut,1)/2,:,m);
    sst_cut(1:size(sst_cut,1)/2,:,m) = sst_cut(size(sst_cut,1)/2+1:end,:,m);
    sst_cut(size(sst_cut,1)/2+1:end,:,m) = ttt;
end

%%
sst_lat = flip(sst_lat);
tt = sst_lon(1:size(sst_cut,1)/2);
sst_lon(1:size(sst_cut,1)/2) = sst_lon(size(sst_cut,1)/2+1:end);
sst_lon(size(sst_cut,1)/2+1:end) = tt+360;
% tttt = ncread(path_index,'time');
sst_t = double(ncread(path_index,'time') + datenum(1870,1,0));
%%
sst_cut(find (sst_cut < -500)) = NaN;
%%

nccreate('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','sst','Dimensions', {'lon',360,'lat',180,'t', 432});
nccreate('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','lon','Dimensions', {'lon',360});
nccreate('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','lat','Dimensions', {'lat',180});
nccreate('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','time','Dimensions', {'t',432});
ncwrite('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','sst', sst_cut);
ncwrite('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','lon', sst_lon);
ncwrite('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','lat', sst_lat);
ncwrite('sst_Amon_FGOALS-f3-L_historical_r1i1p1f1_gr_197901-201412.nc','time', sst_t(1309:1740));

%%

















