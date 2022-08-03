clear all
close all
clc
%%
river_name = 'volga';

%% Surface or full run-off
% mod = "surf";
mod = "full";

save = true;
% save = false;

%%
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman');
global col_arr weights_list list_of_models_hist
col_arr = ['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';...
    '#A2142F';'#00FF00'];
% weights_list = [1, 6, 7];
weights_list = [1, 2, 3, 4, 5, 6, 7, 8];
% figure('Units', 'normalized', 'OuterPosition', [.3 .3 .4 .4]);
%%
% load hist_amur13.07.21.mat% ____________________________________________
% load ssp_amur13.07.21.mat
% load rivers_data_year/hist_amur_29.10.21.mat% ____________________________________________
% load rivers_data_year/ssp_amur_29.10.21.mat
%%
%

load rivers_data_month/pdo(tas)_1979-2014_month_11.02.22.mat
load rivers_data_month/pdo(tas)_2015-2100_month_11.02.22.mat

% load rivers_data_month/nao_cut(psl)_1979-2014_month_18.04.22.mat
load rivers_data_month/nao_cut(psl)_1980-2014_month_22.05.22.mat
load rivers_data_month/nao_cut(psl)_2015-2100_month_18.04.22.mat

% load rivers_data_month/pr_amur_2015-2100_month_25.01.22.mat
% load rivers_data_month/pr_amur_1979-2014_month_25.01.22.mat
% load rivers_data_month/mrro_amur_2015-2100_month_25.01.22.mat
% load rivers_data_month/mrro_amur_1979-2014_month_25.01.22.mat
% load rivers_data_month/mrros_amur_2015-2100_month_25.01.22.mat
% load rivers_data_month/mrros_amur_1979-2014_month_25.01.22.mat

% load rivers_data_month/pr_selenga_2015-2100_month_25.01.22.mat
% load rivers_data_month/pr_selenga_1979-2014_month_25.01.22.mat
% load rivers_data_month/mrro_selenga_2015-2100_month_25.01.22.mat
% load rivers_data_month/mrro_selenga_1979-2014_month_25.01.22.mat
% load rivers_data_month/mrros_selenga_2015-2100_month_25.01.22.mat
% load rivers_data_month/mrros_selenga_1979-2014_month_25.01.22.mat

load rivers_data_month/pr_volga_2015-2100_month_25.01.22.mat
% load rivers_data_month/pr_volga_1979-2014_month_25.01.22.mat
load rivers_data_month/pr_volga_1980-2014_month_22.05.22.mat

% load rivers_data_month/mrro_volga_2015-2100_month_25.01.22.mat
load rivers_data_month/mrro_volga_1980-2014_month_22.05.22.mat
load rivers_data_month/mrro_volga_2015-2100_month_25.01.22.mat

load rivers_data_month/mrros_volga_2015-2100_month_25.01.22.mat
load rivers_data_month/mrros_volga_1980-2014_month_22.05.22.mat

%%
list_of_models_126 = list_of_models_126';
list_of_models_245 = list_of_models_245';
list_of_models_585 = list_of_models_585';
%% just a plot
% figure;
% h = plot(years,R_245,'LineWidth',2);
% legend(list_of_models_245)

% h(1).LineWidth = 2.5;
% все веса считать до 14 года, то есть по данным хисторикал. дельта ка
% (каэф лин тренда и его стандартное отклонение делать по книге
% НумРекФортран стр 686(там другие обозначения). Чуть позже построить корреляционные функции
% сигнала и аппроксимировать экспонентой, найти значение "декремента" -
% когда амплитуда падает в е раз, потом пересчитать дельта к по формуле в
% книге (там будет суммирование не по всему периоду, а по 2*декремент.
%% convert month data to year data to calculate Wm

P_126_hist_y = month_to_year_data(filter_by_model(P_hist_m, list_of_models_126));
R_126_hist_y = month_to_year_data(filter_by_model(R_hist_m, list_of_models_126));
Rs_126_hist_y = month_to_year_data(filter_by_model(Rs_hist_m, list_of_models_126));

P_245_hist_y = month_to_year_data(filter_by_model(P_hist_m, list_of_models_245));
R_245_hist_y = month_to_year_data(filter_by_model(R_hist_m, list_of_models_245));
Rs_245_hist_y = month_to_year_data(filter_by_model(Rs_hist_m, list_of_models_245));

P_585_hist_y = month_to_year_data(filter_by_model(P_hist_m, list_of_models_585));
R_585_hist_y = month_to_year_data(filter_by_model(R_hist_m, list_of_models_585));
Rs_585_hist_y = month_to_year_data(filter_by_model(Rs_hist_m, list_of_models_585));
%%

P_126_ssp_y = month_to_year_data(P_126_21_m);
R_126_ssp_y = month_to_year_data(R_126_21_m);
Rs_126_ssp_y = month_to_year_data(Rs_126_21_m);

P_245_ssp_y = month_to_year_data(P_245_21_m);
R_245_ssp_y = month_to_year_data(R_245_21_m);
Rs_245_ssp_y = month_to_year_data(Rs_245_21_m);

P_585_ssp_y = month_to_year_data(P_585_21_m);
R_585_ssp_y = month_to_year_data(R_585_21_m);
Rs_585_ssp_y = month_to_year_data(Rs_585_21_m);








%% calc Wm
% weights = [w_m,w_tr,w_iav,y_tmp_mean',y_tmp_std',k_tmp',delta_k_tmp',delta_y_tmp_std'];
[P_126_w] = calc_norm_all_models(years_hist, P_126_hist_y );
[R_126_w] = calc_norm_all_models(years_hist, R_126_hist_y);
[Rs_126_w] = calc_norm_all_models(years_hist, Rs_126_hist_y);

[P_245_w] = calc_norm_all_models(years_hist, P_245_hist_y);
[R_245_w] = calc_norm_all_models(years_hist, R_245_hist_y);
[Rs_245_w] = calc_norm_all_models(years_hist, Rs_245_hist_y);

[P_585_w] = calc_norm_all_models(years_hist, P_585_hist_y);
[R_585_w] = calc_norm_all_models(years_hist, R_585_hist_y);
[Rs_585_w] = calc_norm_all_models(years_hist, Rs_585_hist_y);

%% Calculating W_nv
% P_w_nv = calc_w_nv(P_hist_m, pdo_hist);
% R_w_nv = calc_w_nv(R_hist_m, pdo_hist);
% Rs_w_nv = calc_w_nv(Rs_hist_m, pdo_hist);

P_w_nv = calc_w_nv(P_hist_m, nao_hist_cut);
R_w_nv = calc_w_nv(R_hist_m, nao_hist_cut);
Rs_w_nv = calc_w_nv(Rs_hist_m, nao_hist_cut);

weights_nv = P_w_nv.*R_w_nv/sum(P_w_nv.*R_w_nv);
weights_nv_s = P_w_nv.*Rs_w_nv/sum(P_w_nv.*Rs_w_nv);

%% Calculating Y_f

P_hist_f = calc_y_f(P_hist_m, nao_hist_cut);
R_hist_f = calc_y_f(R_hist_m, nao_hist_cut);
Rs_hist_f = calc_y_f(Rs_hist_m, nao_hist_cut);

%%
% P_126_ssp_f = calc_y_f(P_126_21_m, pdo_126_21);
% R_126_ssp_f = calc_y_f(R_126_21_m, pdo_126_21);
% Rs_126_ssp_f = calc_y_f(Rs_126_21_m, pdo_126_21);
% 
% P_245_ssp_f = calc_y_f(P_245_21_m, pdo_245_21);
% R_245_ssp_f = calc_y_f(R_245_21_m, pdo_245_21);
% Rs_245_ssp_f = calc_y_f(Rs_245_21_m, pdo_245_21);
% 
% P_585_ssp_f = calc_y_f(P_585_21_m, pdo_585_21);
% R_585_ssp_f = calc_y_f(R_585_21_m, pdo_585_21);
% Rs_585_ssp_f = calc_y_f(Rs_585_21_m, pdo_585_21);

P_126_ssp_f = calc_y_f(P_126_21_m, nao_126_21_cut);
R_126_ssp_f = calc_y_f(R_126_21_m, nao_126_21_cut);
Rs_126_ssp_f = calc_y_f(Rs_126_21_m, nao_126_21_cut);

P_245_ssp_f = calc_y_f(P_245_21_m, nao_245_21_cut);
R_245_ssp_f = calc_y_f(R_245_21_m, nao_245_21_cut);
Rs_245_ssp_f = calc_y_f(Rs_245_21_m, nao_245_21_cut);

P_585_ssp_f = calc_y_f(P_585_21_m, nao_585_21_cut);
R_585_ssp_f = calc_y_f(R_585_21_m, nao_585_21_cut);
Rs_585_ssp_f = calc_y_f(Rs_585_21_m, nao_585_21_cut);

%%
%%
%% iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
%%






P_126_hist_f = month_to_year_data(filter_by_model(P_hist_f, list_of_models_126));
R_126_hist_f = month_to_year_data(filter_by_model(R_hist_f, list_of_models_126));
Rs_126_hist_f = month_to_year_data(filter_by_model(Rs_hist_f, list_of_models_126));

P_245_hist_f = month_to_year_data(filter_by_model(P_hist_f, list_of_models_245));
R_245_hist_f = month_to_year_data(filter_by_model(R_hist_f, list_of_models_245));
Rs_245_hist_f = month_to_year_data(filter_by_model(Rs_hist_f, list_of_models_245));

P_585_hist_f = month_to_year_data(filter_by_model(P_hist_f, list_of_models_585));
R_585_hist_f = month_to_year_data(filter_by_model(R_hist_f, list_of_models_585));
Rs_585_hist_f = month_to_year_data(filter_by_model(Rs_hist_f, list_of_models_585));
%%

P_126_ssp_y_f = month_to_year_data(P_126_ssp_f);
R_126_ssp_y_f = month_to_year_data(R_126_ssp_f);
Rs_126_ssp_y_f = month_to_year_data(Rs_126_ssp_f);

P_245_ssp_y_f = month_to_year_data(P_245_ssp_f);
R_245_ssp_y_f = month_to_year_data(R_245_ssp_f);
Rs_245_ssp_y_f = month_to_year_data(Rs_245_ssp_f);

P_585_ssp_y_f = month_to_year_data(P_585_ssp_f);
R_585_ssp_y_f = month_to_year_data(R_585_ssp_f);
Rs_585_ssp_y_f = month_to_year_data(Rs_585_ssp_f);








%% calc Wm with Y_f
% weights = [w_m,w_tr,w_iav,y_tmp_mean',y_tmp_std',k_tmp',delta_k_tmp',delta_y_tmp_std'];
[P_126_w_f] = calc_norm_all_models(years_hist, P_126_hist_f );
[R_126_w_f] = calc_norm_all_models(years_hist, R_126_hist_f);
[Rs_126_w_f] = calc_norm_all_models(years_hist, Rs_126_hist_f);

[P_245_w_f] = calc_norm_all_models(years_hist, P_245_hist_f);
[R_245_w_f] = calc_norm_all_models(years_hist, R_245_hist_f);
[Rs_245_w_f] = calc_norm_all_models(years_hist, Rs_245_hist_f);

[P_585_w_f] = calc_norm_all_models(years_hist, P_585_hist_f);
[R_585_w_f] = calc_norm_all_models(years_hist, R_585_hist_f);
[Rs_585_w_f] = calc_norm_all_models(years_hist, Rs_585_hist_f);

%%
P_126_w_f(2:end,1) = P_126_w(2:end,1);
R_126_w_f(2:end,1) = R_126_w(2:end,1);
Rs_126_w_f(2:end,1) = Rs_126_w(2:end,1);

P_245_w_f(2:end,1) = P_245_w(2:end,1);
R_245_w_f(2:end,1) = R_245_w(2:end,1);
Rs_245_w_f(2:end,1) = Rs_245_w(2:end,1);

P_585_w_f(2:end,1) = P_585_w(2:end,1);
R_585_w_f(2:end,1) = R_585_w(2:end,1);
Rs_585_w_f(2:end,1) = Rs_585_w(2:end,1);

%% Now weights = [w_m,w_tr,w_iav,w_nv,y_tmp_std',k_tmp',delta_k_tmp',delta_y_tmp_std'];

P_126_w_f(:,4) = filter_by_model(P_w_nv', list_of_models_126);
R_126_w_f(:,4) = filter_by_model(R_w_nv', list_of_models_126);
Rs_126_w_f(:,4) = filter_by_model(Rs_w_nv', list_of_models_126);

P_245_w_f(:,4) = filter_by_model(P_w_nv', list_of_models_245);
R_245_w_f(:,4) = filter_by_model(R_w_nv', list_of_models_245);
Rs_245_w_f(:,4) = filter_by_model(Rs_w_nv', list_of_models_245);

P_585_w_f(:,4) = filter_by_model(P_w_nv', list_of_models_585);
R_585_w_f(:,4) = filter_by_model(R_w_nv', list_of_models_585);
Rs_585_w_f(:,4) = filter_by_model(Rs_w_nv', list_of_models_585);

%%

% weights = [w_m, w_tr, w_iav, w_nv, w_r, w_p, w_all];
weights_126 = calc_weights_f(P_126_w_f,R_126_w_f);
weights_126_s = calc_weights_f(P_126_w_f,Rs_126_w_f);

weights_245 = calc_weights_f(P_245_w_f,R_245_w_f);
weights_245_s = calc_weights_f(P_245_w_f,Rs_245_w_f);

weights_585 = calc_weights_f(P_585_w_f,R_585_w_f);
weights_585_s = calc_weights_f(P_585_w_f,Rs_585_w_f);















%%

% figure;
% plot(P_hist_m(1,:));
% hold on
% plot(pdo_hist(1,:)*k_tmp);


%%
% z = 1979:2014;
% plot(z, output2,'o')
% hold on
% plot(z,P_585(2,1:end-1));
% plot(z,month_to_year_data(P_hist(1, :)),'*');
% 
% %%
% a = month_to_year_data(P_hist);
% plot(z, a)
% figure;
% 
% plot(z, P_585(:,1:end-1))
%%

% plot(years_hist, P_126(:,1:end-1))
% figure;
% 
% plot(years_hist, P_126_hist_y)
%%













%%



%%
% weights = [w_m,w_tr,w_iav,w_r,w_p,w_all];
% weights_126 = calc_weights(P_126_w,R_126_w);
% weights_126_s = calc_weights(P_126_w,Rs_126_w);
% 
% weights_245 = calc_weights(P_245_w,R_245_w);
% weights_245_s = calc_weights(P_245_w,Rs_245_w);
% 
% weights_585 = calc_weights(P_585_w,R_585_w);
% weights_585_s = calc_weights(P_585_w,Rs_585_w);
%% weights for simple mean 
weights_126_s(:,end+1) = 1/(size(weights_126_s,1)-1);
weights_245_s(:,end+1) = 1/(size(weights_245_s,1)-1);
weights_585_s(:,end+1) = 1/(size(weights_585_s,1)-1);

weights_126(:,end+1) = 1/(size(weights_126,1)-1);
weights_245(:,end+1) = 1/(size(weights_245,1)-1);
weights_585(:,end+1) = 1/(size(weights_585,1)-1);

%%
global list_of_weights 
list_of_weights = {'W_m','W_t_r','W_i_a_v','W_n_v','W_r','W_p','W_a_l_l','W_0'};
list_of_scenarios = {'SSP1-2.6','SSP2-4.5','SSP5-8.5'};

%%

% P_126(:,end+1) = P_126_tmp;
% R_126(:,end+1) = R_126_tmp;
% Rs_126(:,end+1) = Rs_126_tmp;
% 
% P_245(:,end+1) = P_245_tmp;
% R_245(:,end+1) = R_245_tmp;
% Rs_245(:,end+1) = Rs_245_tmp;
% 
% P_585(:,end+1) = P_585_tmp;
% R_585(:,end+1) = R_585_tmp;
% Rs_585(:,end+1) = Rs_585_tmp;
% 
% years(end+1) = 2015;

%%
% clc
if mod == "full"
    [means_126,u_z_126] = calc_means(weights_126,R_126_hist_f(2:end,:));
    [means_245,u_z_245] = calc_means(weights_245,R_245_hist_f(2:end,:));
    [means_585,u_z_585] = calc_means(weights_585,R_585_hist_f(2:end,:));

    [means_126_21,u_z_126_21] = calc_means(weights_126,R_126_ssp_y_f);
    [means_245_21,u_z_245_21] = calc_means(weights_245,R_245_ssp_y_f);
    [means_585_21,u_z_585_21] = calc_means(weights_585,R_585_ssp_y_f);
elseif mod == "surf"
    [means_126,u_z_126] = calc_means(weights_126_s,Rs_hist_f(2:end,:));
    [means_245,u_z_245] = calc_means(weights_245_s,Rs_hist_f(2:end,:));
    [means_585,u_z_585] = calc_means(weights_585_s,Rs_hist_f(2:end,:));

    [means_126_21,u_z_126_21] = calc_means(weights_126_s,Rs_126_ssp_y_f);
    [means_245_21,u_z_245_21] = calc_means(weights_245_s,Rs_245_ssp_y_f);
    [means_585_21,u_z_585_21] = calc_means(weights_585_s,Rs_585_ssp_y_f);
else
    disp('mod error')
end


%% calc std
period = 11; % ----------------------------------------------------------------------------------------------------------------------------

years_full = [years_hist,years_ssp];
yrs = years_full(1) + period/2 :1: years_full(end) - period/2;

std_126_s = calc_std(Rs_126_hist_y,Rs_126_ssp_y,period);
std_245_s = calc_std(Rs_245_hist_y,Rs_245_ssp_y,period);
std_585_s = calc_std(Rs_585_hist_y,Rs_585_ssp_y,period);

std_126 = calc_std(R_126_hist_y,R_126_ssp_y,period);
std_245 = calc_std(R_245_hist_y,R_245_ssp_y,period);
std_585 = calc_std(R_585_hist_y,R_585_ssp_y,period);


[std_126_s_means,uuu1] = calc_means(weights_126_s,std_126_s);
[std_245_s_means,uuu2] = calc_means(weights_245_s,std_245_s);
[std_585_s_means,uuu3] = calc_means(weights_585_s,std_585_s);
[std_126_means,uuu1] = calc_means(weights_126,std_126);
[std_245_means,uuu2] = calc_means(weights_245,std_245);
[std_585_means,uuu3] = calc_means(weights_585,std_585);
%% plot std prep
yrs_etalon = years_hist(1) + period/2 :1: years_hist(end) - period/2;
var = Rs_126_hist_y(1,:);
len_of_std = numel(var);
std_out = zeros(1,len_of_std-period);

for kk = 1:len_of_std - period
    tmp = var(kk : kk + period - 1);
%         std_out(model_n,kk) = std(tmp,0,2);
    std_etalon_s(kk) = std(tmp);
end


%% plot std
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');
if mod == "surf"
    var_tmp(1,:,:) = std_126_s_means;
    var_tmp(2,:,:) = std_245_s_means;
    var_tmp(3,:,:) = std_585_s_means;
elseif mod == "full"
    var_tmp(1,:,:) = std_126_means;
    var_tmp(2,:,:) = std_245_means;
    var_tmp(3,:,:) = std_585_means;
else
    display("std error");
end
    
for pp = 1:3
    figure('Units', 'normalized', 'OuterPosition', [.2 .2 .28 .4]);%win
%     figure('Units', 'normalized', 'OuterPosition', [.2 .2 .5 .6]);%mac
    hold on;
    data_for_plot(:,:) = var_tmp(pp,:,:);
    for mm = 1:8
        plot(yrs,data_for_plot(mm,:),'LineWidth',2,'Color',col_arr(mm,:));
    end
    xlabel('год');
    ylabel('км^3/год');
    grid on;
    plot(yrs_etalon,std_etalon_s,'LineWidth',3,'Color','k');
    title(list_of_scenarios(pp));
%     lgd = legend(p1,list_of_weights);
%     legend(list_of_weights);
    axx = gca;
    axx.XLim = [1977 2102];
    if save == true
        exportgraphics(axx,[string(river_name)+"_"+string(list_of_scenarios(pp))+"_std_"+string(period)+"_years_"+mod+".png"])
    end
end
% close all;
% clear var_tmp data_for_plot std_etalon_s
%% Plot Weights prep
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');

w_585_marker = ones(1,size(weights_585_s,1)-1);
for tt = 1:size(weights_585_s,1) - 1
    if any(string(list_of_models_585(tt)) == string(list_of_models_126))
        w_126_marker(tt) = 1;
    else
        w_126_marker(tt) = NaN;
    end
    
    if any(string(list_of_models_585(tt)) == string(list_of_models_245))
        w_245_marker(tt) = 1;
    else
        w_245_marker(tt) = NaN;
    end    
end





%% Plot Weights

if mod == "surf"
    ww = weights_585_s;
elseif mod == "full"
    ww = weights_585;
else
    display("weights error");
end



for nn = 1:7
    figure('Units', 'normalized', 'OuterPosition', [.2 .2 .11 .26]);%win
%     figure('Units', 'normalized', 'OuterPosition', [.2 .2 .18 .38]);%mac
    w_num = nn;
%     X = categorical(list_of_models_585);
%     X = reordercats(X,list_of_models_585);
   
%     bar(ww(2:end,w_num),[0.5 0.5 0.5]);
    b = bar(ww(2:end,w_num));   
    b.FaceColor = [0.5 0.5 0.5];
    coef = max(ww(2:end,w_num));
    
%     bar(X,ww(2:end,w_num));
    % set(gca, 'YScale', 'log');
    % xlabel('Model');
    title(list_of_weights(w_num));

    grid on;
    xlim = get(gca,'xlim');
    hold on
    size_w = size(ww,1);

    plot(xlim,[1/(size_w-1) 1/(size_w-1)],'LineWidth',2,'Color','#4DBEEE')
%     plot(xlim,[1/(15) 1/(15)],'LineWidth',2,'Color','#77AC30')    
%     plot(xlim,[1/(17) 1/(17)],'LineWidth',2,'Color','#A2142F')
    
    ax = gca;
    ax.FontSmoothing = 'on';
%     BX = get(gca,'XTick');
    ax.XTick = [3,6,9,12,15,18];
    plot(w_585_marker*1.1*coef,'o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE');
    plot(w_245_marker*1.2*coef,'o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30');    
    plot(w_126_marker*1.3*coef,'o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F');      
    
    ax.YLim = [0 coef*1.4];
    
%     
%     BY = get(gca,'YTick');
%     xlabel('lol','Position',[BX(size(BX,2)) BY(1)]);
%     ylabel(list_of_weights(w_num),'Rotation',0,'Position',[BX(1) BY(size(BY,2))]);
%     text(1,1,'$\varphi$','FontSize',14,'Interpreter', 'latex')
    if save == true
        exportgraphics(ax, [string(river_name)+"_"+string(list_of_weights(w_num))+"_"+mod+".png"])
    end
   
end
% close all;

%%


% set(0,'DefaultAxesFontSize',24,'DefaultAxesFontName','Times New Roman');
% set(0,'DefaultTextFontSize',24,'DefaultTextFontName','Times New Roman');
% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .6 .6]);
% % figure;
% hold on;
% % for ii = 1:6
% %     h = plot(years_21,means_126(ii,:),'LineWidth',2);
% %     
% % %     errorbar(years,means(ii,:),u_z(ii,:));
% % end 
% plot([years,years_21],[means_585,means_585_21],'LineWidth',2);
% plot(years,Rs_585(1,:),'LineWidth',3,'Color','k');
% plot([years,years_21],[Rs_585(2:end,:),Rs_585_21],'Color',[0.5 0.5 0.5]);
% 
% 
% 
% 
% legend([list_of_weights,{'observ.'}]);
% % plot(years,means);
% title('SSP 585');
% xlabel('years');
% ylabel('Total runoff (mrros), km^3/year');
%% Plot means and std
if mod == "surf"
    Plot_res(years_hist,years_ssp,Rs_585_hist_y,Rs_585_ssp_y,means_585,means_585_21,u_z_585,u_z_585_21);
%     title('SSP5-8.5');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_585_"+mod+".png"]);
    end
    
    Plot_res(years_hist,years_ssp,Rs_245_hist_y,Rs_245_ssp_y,means_245,means_245_21,u_z_245,u_z_245_21);
%     title('SSP2-4.5');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_245_"+mod+".png"]);
    end
    

    Plot_res(years_hist,years_ssp,Rs_126_hist_y,Rs_126_ssp_y,means_126,means_126_21,u_z_126,u_z_126_21);
%     title('SSP1-2.6');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_126_"+mod+".png"]);
    end
    
%     Plot_res_hist(years,years_21,Rs_585,Rs_585_21,means_585,means_585_21,u_z_585,u_z_585_21);
%     title('Hist');
%     ax = gca;
%     if save == true
%         exportgraphics(ax,[string(river_name)+"_mean_hist_"+mod+".png"]);
%     end
    
elseif  mod == "full"
    Plot_res(years_hist,years_ssp,R_585_hist_y,R_585_ssp_y,means_585,means_585_21,u_z_585,u_z_585_21);
%     title('SSP5-8.5');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_585_"+mod+".png"]);
    end
    

    Plot_res(years_hist,years_ssp,R_245_hist_y,R_245_ssp_y,means_245,means_245_21,u_z_245,u_z_245_21);
%     title('SSP2-4.5');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_245_"+mod+".png"]);
    end
    

    Plot_res(years_hist,years_ssp,R_126_hist_y,R_126_ssp_y,means_126,means_126_21,u_z_126,u_z_126_21);
%     title('SSP1-2.6');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_126_"+mod+".png"]);
    end
    
%     Plot_res_hist(years,years_21,R_585,R_585_21,means_585,means_585_21,u_z_585,u_z_585_21);
%     title('Hist');
%     ax = gca;
%     if save == true
%          exportgraphics(ax,[string(river_name)+"_mean_hist_"+mod+".png"]);
%     end
end


% 
% close all;
%% save data
% amur_126_means = [means_126,means_126_21];
% amur_245_means = [means_245,means_245_21];
% amur_585_means = [means_585,means_585_21];
% 
% amur_years = [years,years_21];
% %%
% 
% save('amur.mat','amur_126_means','amur_245_means','amur_585_means','amur_years');

%% Plot smooth means and std
% 
% save = true;
if mod == "surf"
    Plot_smooth(years_hist,years_ssp,means_585,means_585_21,u_z_585,u_z_585_21,0.2);
%     title('SSP5-8.5');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_585_"+mod+"_smooth.png"]);
    end
    
    Plot_smooth(years_hist,years_ssp,means_245,means_245_21,u_z_245,u_z_245_21,0.2);
%     title('SSP2-4.5');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_245_"+mod+"_smooth.png"]);
    end
    

    Plot_smooth(years_hist,years_ssp,means_126,means_126_21,u_z_126,u_z_126_21,0.2);
%     title('SSP1-2.6');
    ax = gca;
    if save == true
        exportgraphics(ax,[string(river_name)+"_mean_126_"+mod+"_smooth.png"]);
    end

elseif  mod == "full"
    
    Plot_smooth(years_hist,years_ssp,means_585,means_585_21,u_z_585,u_z_585_21,0.2);
%     title('SSP5-8.5');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_585_"+mod+"_smooth.png"]);
    end
    

    Plot_smooth(years_hist,years_ssp,means_245,means_245_21,u_z_245,u_z_245_21,0.2);
%     title('SSP2-4.5');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_245_"+mod+"_smooth.png"]);
    end
    

    Plot_smooth(years_hist,years_ssp,means_126,means_126_21,u_z_126,u_z_126_21,0.2);
%     title('SSP1-2.6');
    ax = gca;
    if save == true
         exportgraphics(ax,[string(river_name)+"_mean_126_"+mod+"_smooth.png"]);
    end

end

% close all;
%%
% years_full = [years,years_21];
% for bb = 1:6
%     [mean_126(bb),std_126(bb),k_126(bb),delta_k_126(bb),delta_y_126(bb)] = lin_reg(years_full,[means_126(bb,:),means_126_21(bb,:)]);
%     [mean_245(bb),std_245(bb),k_245(bb),delta_k_245(bb),delta_y_245(bb)] = lin_reg(years_full,[means_245(bb,:),means_245_21(bb,:)]);    
%     [mean_585(bb),std_585(bb),k_585(bb),delta_k_585(bb),delta_y_585(bb)] = lin_reg(years_full,[means_585(bb,:),means_585_21(bb,:)]);    
% end




%%
% set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
% set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');
% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .4 .6]);
% plot([years,years_21],[u_z_126,u_z_126_21],'LineWidth',2);
% grid on;
% xlabel('год');
% ylabel('км^3/год');
% legend(list_of_weights);


%%
% figure
% plot(normpdf(1:w_length,P_585_m,P_585_std))
% a = normpdf(1:w_length,P_585_m,P_585_std);

%%

% tmp = P_585(1,:);
% cor_tmp = ifft(fft(tmp).*conj(fft(tmp)));
% PlotCorr(cor_tmp,1);
% 
% figure;
% 
% plot();
% a=xcorr(tmp,tmp,20,'coeff');

%%

% scatter(P_hist_m(1,:), nao_hist_cut(1,:));
% %%
plot(P_hist_m(1,:)/max(P_hist_m(1,:)));
hold on
plot(pdo_hist(1,:)/max(pdo_hist(1,:)));
%%
plot(R_hist_m(1,:)/max(R_hist_m(1,:)));
hold on
plot(pdo_hist(1,:)/max(pdo_hist(1,:)));
%%
[a,b] = my_lin_reg(nao_hist_cut(1,:), P_hist_m(1,:));
ttt = linspace(1980, 2014, 420);

plot(ttt, nao_hist_cut);
hold on
plot(ttt, nao_hist_cut(1,:),'LineWidth',3,'Color','k')
figure;
plot(ttt, P_hist_m);
%% std calc
len_tmp = 344;
for i = 1:3
    for j = 1 : size(P_126_ssp_f(1,:))
        P_std_orig_126(j,i) = std(P_126_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_orig_126(j,i) = std(R_126_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        
        P_std_orig_245(j,i) = std(P_245_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_orig_245(j,i) = std(R_245_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        
        P_std_orig_585(j,i) = std(P_585_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_orig_585(j,i) = std(R_585_ssp_f(j,(i-1)*len_tmp+1:i*len_tmp));
        
        P_std_eof_126(j,i) = std(P_126_21_m(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_eof_126(j,i) = std(R_126_21_m(j,(i-1)*len_tmp+1:i*len_tmp));
        
        P_std_eof_245(j,i) = std(P_245_21_m(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_eof_245(j,i) = std(R_245_21_m(j,(i-1)*len_tmp+1:i*len_tmp));
        
        P_std_eof_585(j,i) = std(P_585_21_m(j,(i-1)*len_tmp+1:i*len_tmp));
        R_std_eof_585(j,i) = std(R_585_21_m(j,(i-1)*len_tmp+1:i*len_tmp));    
        
    end
    
end
%% relative std

P_rel_std_126 = (1 - P_std_orig_126./P_std_eof_126)*100;
R_rel_std_126 = (1 - R_std_orig_126./R_std_eof_126)*100;

P_rel_std_245 = (1 - P_std_orig_245./P_std_eof_245)*100;
R_rel_std_245 = (1 - R_std_orig_245./R_std_eof_245)*100;

P_rel_std_585 = (1 - P_std_orig_585./P_std_eof_585)*100;
R_rel_std_585 = (1 - R_std_orig_585./R_std_eof_585)*100;

%%
decs = 1:3;
hold on
plot(decs, P_rel_std_126,'LineWidth',3);
plot(decs, R_rel_std_126,'LineWidth',3);
plot(decs, P_rel_std_245,'LineWidth',3);
plot(decs, R_rel_std_245,'LineWidth',3);
plot(decs, P_rel_std_585,'LineWidth',3);
plot(decs, R_rel_std_585,'LineWidth',3);


%%
% 
% figure;
% plot_data_and_trend(P_585(1,:));

% count = 1;
% data = P_585;
% yrs = years;
function [std_out] = calc_std(var_hist,var_ssp,period)
n_of_models = size(var_hist,1)-1;
var = [var_hist(2:end,:),var_ssp];
len_of_std = numel(var(1,:));
std_out = zeros(n_of_models,len_of_std-period);

for model_n = 1:n_of_models
    for kk = 1:len_of_std - period
        tmp = var(model_n,kk : kk + period - 1);
%         std_out(model_n,kk) = std(tmp,0,2);
        std_out(model_n,kk) = std(tmp);
    end
end

end

function Plot_res(years_hist,years_ssp,var_hist,var_ssp,mean_hist,mean_ssp,u_hist,u_ssp)
global list_of_weights col_arr weights_list
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');


% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .3 .5]);
% hold on;
% for mm = 1:7
%     plot([years_hist,years_ssp],[u_hist(mm,:),u_ssp(mm,:)],'LineWidth',2,'Color',col_arr(mm,:));
% end
% xlabel('год');
% ylabel('км^3/год');
% grid on;

figure('Units', 'normalized', 'OuterPosition', [.2 .2 .28 .4]);%win
% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .5 .6]);%mac
hold on;
p1 = plot([years_hist,years_ssp],[var_hist(2:end,:),var_ssp],'Color',[0.5 0.5 0.5]);
% p1 = plot(years_ssp,var_ssp,'Color',[0.5 0.5 0.5]);

for nn = weights_list
    p1 = plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:));
%     p1 = plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:),'DisplayName',char(list_of_weights(nn)));
%     p1 = plot(years_ssp, mean_ssp(nn,:),'LineWidth',2,'Color',col_arr(nn,:),'DisplayName',char(list_of_weights(nn)));
end

% plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k','DisplayName','эталон.');
plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k');

% lgd = legend('boxoff');
% lgd.FontSize = 20;
% lgd = legend('Orientation','horizontal');
% plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k');

% xlabel('years');
% ylabel('surface runoff (mrros), km^3/year');
xlabel('год');
ylabel('км^3/год');
grid on;
ax = gca;
ax.FontSmoothing = 'on';
ax.XLim = [years_hist(1)-2 years_ssp(end)+2];
% ax.XLim = [years_ssp(1)-1 years_ssp(end)+1];
% ax.YLim = [0 90];
end


function Plot_res_hist(years_hist,years_ssp,var_hist,var_ssp,mean_hist,mean_ssp,u_hist,u_ssp)
global list_of_weights col_arr weights_list
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');


% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .3 .5]);
% hold on;
% for mm = 1:7
%     plot([years_hist,years_ssp],[u_hist(mm,:),u_ssp(mm,:)],'LineWidth',2,'Color',col_arr(mm,:));
% end
% xlabel('год');
% ylabel('км^3/год');
% grid on;

figure('Units', 'normalized', 'OuterPosition', [.2 .2 .28 .4]);
hold on;
% p1 = plot([years_hist,years_ssp],[var_hist(2:end,:),var_ssp],'Color',[0.5 0.5 0.5]);
p1 = plot(years_hist,var_hist(2:end,:),'Color',[0.5 0.5 0.5]);

for nn = weights_list
%     p1 = plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:));
%     p1 = plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:),'DisplayName',char(list_of_weights(nn)));
    p1 = plot(years_hist, mean_hist(nn,:),'LineWidth',2,'Color',col_arr(nn,:),'DisplayName',char(list_of_weights(nn)));
end

% plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k','DisplayName','эталон.');
plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k');

% lgd = legend('boxoff');
% lgd.FontSize = 20;
% lgd = legend('Orientation','horizontal');
% plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k');

% xlabel('years');
% ylabel('surface runoff (mrros), km^3/year');
xlabel('год');
ylabel('км^3/год');
grid on;
ax = gca;
ax.FontSmoothing = 'on';
ax.XLim = [years_hist(1)-1 years_hist(end)+1];
% ax.XLim = [years_ssp(1)-1 years_ssp(end)+1];
% ax.YLim = [0 90];
end


function Plot_smooth(years_hist,years_ssp,mean_hist,mean_ssp,u_hist,u_ssp,sm_factor)
global list_of_weights col_arr
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman');


% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .4 .6]);
% hold on;
% for mm = 1:6
%     
% %     plot([years_hist,years_ssp],[u_hist(mm,:),u_ssp(mm,:)],'LineWidth',2,'Color',col_arr(mm,:));
%     plot([years_hist';years_ssp'],smooth([years_hist,years_ssp],[u_hist(mm,:),u_ssp(mm,:)],sm_factor),'LineWidth',2,'Color',col_arr(mm,:));
% end
% xlabel('год');
% ylabel('км^3/год');
% grid on;

figure('Units', 'normalized', 'OuterPosition', [.2 .2 .28 .4]);%win
% figure('Units', 'normalized', 'OuterPosition', [.2 .2 .5 .6]);%mac

hold on;

for nn = 1:6
%     plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:));
    plot([years_hist';years_ssp'],smooth([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],sm_factor),'LineWidth',2,'Color',col_arr(nn,:));
%     p1 = plot([years_hist,years_ssp],[mean_hist(nn,:),mean_ssp(nn,:)],'LineWidth',2,'Color',col_arr(nn,:),'DisplayName',char(list_of_weights(nn)));
end

% lgd = legend(p1,list_of_weights);
% lgd = legend('boxoff');
% lgd.FontSize = 16;
% lgd = legend('Orientation','horizontal');
% plot(years_hist,var_hist(1,:),'LineWidth',3,'Color','k');

% xlabel('years');
% ylabel('surface runoff (mrros), km^3/year');
xlabel('год');
ylabel('км^3/год');
grid on;
ax = gca;
ax.FontSmoothing = 'on';
end



function [means,u_z] = calc_means(weights,vars)

    number_of_models = size(vars,1);
    n_of_weights = size(weights,2);
    means = zeros(n_of_weights,size(vars,2));
    u_z = means;
    % calc mean
    for z_count = 1 : n_of_weights

        for model_n = 1 : number_of_models
            mean_tmp = vars(model_n,:)*weights(model_n + 1,z_count);
            means(z_count,:) = means(z_count,:) + mean_tmp;
        end
    %     plot(means(z_count,:))
    end

    % calc std
    for z_count = 1 : n_of_weights
        std_tmp = zeros(1,size(vars,2));
        for model_n = 1 : number_of_models
            var_tmp = vars(model_n,:);
            sigma_tmp = std(var_tmp,0,2);
            std_tmp = std_tmp + (var_tmp.^2 + sigma_tmp^2) * weights(model_n + 1,z_count);

        end
       u_z(z_count,:) = sqrt(std_tmp - means(z_count,:).^2);
    end
end


function [weights] = calc_weights(P_w,R_w)
w_m = P_w(:,1).*R_w(:,1);
w_tr = P_w(:,2).*R_w(:,2);
w_iav = P_w(:,3).*R_w(:,3);
w_r = R_w(:,1).*R_w(:,2).*R_w(:,3);
w_p = P_w(:,1).*P_w(:,2).*P_w(:,3);
w_all = w_r.*w_p;

w_m = w_m./sum(w_m);
w_tr = w_tr./sum(w_tr);
w_iav = w_iav./sum(w_iav);
w_r = w_r./sum(w_r);
w_p = w_p./sum(w_p);
w_all = w_all./sum(w_all);

weights = [w_m,w_tr,w_iav,w_r,w_p,w_all];
end

function [weights] = calc_weights_f(P_w,R_w)
w_m = P_w(:,1).*R_w(:,1);
w_tr = P_w(:,2).*R_w(:,2);
w_iav = P_w(:,3).*R_w(:,3);
w_nv = P_w(:,4).*R_w(:,4);

w_r = R_w(:,1).*R_w(:,2).*R_w(:,3).*R_w(:,4);
w_p = P_w(:,1).*P_w(:,2).*P_w(:,3).*P_w(:,4);
w_all = w_r.*w_p.*w_nv;

w_m = w_m./sum(w_m);
w_tr = w_tr./sum(w_tr);
w_iav = w_iav./sum(w_iav);
w_nv = w_nv./sum(w_nv);

w_r = w_r./sum(w_r);
w_p = w_p./sum(w_p);
w_all = w_all./sum(w_all);

weights = [w_m, w_tr, w_iav, w_nv, w_r, w_p, w_all];
end




function [weights] = calc_norm_all_models(yrs, data)
num_of_models = size(data,1);

w_m = zeros(num_of_models,1);
w_tr = zeros(num_of_models,1);
w_iav = zeros(num_of_models,1);
teta = (2/(size(data,2)-1))^(1/4);
[y_tmp_mean(1),y_tmp_std(1),k_tmp(1),delta_k_tmp(1),delta_y_tmp_std(1)] = lin_reg(yrs,data(1,:));
w_m(1) = [];
w_tr(1) = [];
w_iav(1) = [];




    for count = 2 : num_of_models
        [y_tmp_mean(count),y_tmp_std(count),k_tmp(count),delta_k_tmp(count),delta_y_tmp_std(count)] = lin_reg(yrs,data(count,:));

        w_m(count) = exp(-((y_tmp_mean(count)-y_tmp_mean(1))^2)/(2*y_tmp_std(1)^2));
        w_tr(count) = exp(-((k_tmp(count)-k_tmp(1))^2)/(2*delta_k_tmp(1)^2));
        w_iav(count) = exp(-((delta_y_tmp_std(count)-delta_y_tmp_std(1))^2)/(2*(teta*delta_y_tmp_std(count))^2));
    end
    weights = [w_m,w_tr,w_iav,y_tmp_mean',y_tmp_std',k_tmp',delta_k_tmp',delta_y_tmp_std'];
end


function [y_mean,y_std,k,delta_k,delta_y_std] = lin_reg(x,y)
%           Y=A+B*X
%regcoef - коэффициент регрессии (который у нас k),
%regcoefstd - стандартное отклонение для regcoef (у нас оно delta k),
% regr - коэффициент корреляции для регрессии (он в Вашей схеме не нужен),
% regss - статистическая значимость (тоже сейчас не нужно).

x_mean = mean(x);
y_mean = mean(y);
y_std = std(y,0,2);
x_variance = var(x);
y_variance = var(y);
xy_mean = mean(x.*y);
regcoef = (xy_mean-x_mean*y_mean)/x_variance;
regr = regcoef*sqrt(x_variance/y_variance);
rrr = xcorr(x,x,1,'coeff');
rx = rrr(1);
rrr = xcorr(y,y,1,'coeff');

ry = rrr(1);

cor_tmpx = ifft(fft(x).*conj(fft(x)));
cor_tmpy = ifft(fft(y).*conj(fft(y)));

coefdof = (1-abs(rx*ry))/(1+abs(rx*ry));
eqn = length(x)*coefdof;
regcoefstd = sqrt((y_variance/x_variance-regcoef^2)/eqn);
dof = (length(x)-2)*coefdof;
regss = 0;
if(isnan(regcoefstd))
    disp('k_std NaN!');
    [regcoefstd x_mean y_mean x_variance y_variance]
end
  
k = regcoef;
delta_k = regcoefstd;
delta_y = detrend(y);
delta_y_std = std(delta_y,0,2);
  
end

function [k,delta_k] = my_lin_reg(x,y)
%           Y=A+B*X
%regcoef - коэффициент регрессии (который у нас k),
%regcoefstd - стандартное отклонение для regcoef (у нас оно delta k),
% regr - коэффициент корреляции для регрессии (он в Вашей схеме не нужен),
% regss - статистическая значимость (тоже сейчас не нужно).

x_mean = mean(x);
y_mean = mean(y);
y_std = std(y,0,2);
x_variance = var(x);
y_variance = var(y);
xy_mean = mean(x.*y);
regcoef = (xy_mean-x_mean*y_mean)/(mean(x.^2)-x_mean.^2);
regr = regcoef*sqrt(x_variance/y_variance);
rrr = xcorr(x,x,1,'coeff');
rx = rrr(1);
rrr = xcorr(y,y,1,'coeff');

ry = rrr(1);

cor_tmpx = ifft(fft(x).*conj(fft(x)));
cor_tmpy = ifft(fft(y).*conj(fft(y)));

coefdof = (1-abs(rx*ry))/(1+abs(rx*ry));
eqn = length(x)*coefdof;
regcoefstd = sqrt((y_variance/x_variance-regcoef^2)/eqn);

if(isnan(regcoefstd))
    disp('k_std NaN!');
    [regcoefstd x_mean y_mean x_variance y_variance]
end
  
k = regcoef;
delta_k = regcoefstd;

end

function[output] = month_to_year_data(mon_data)
    if size(mon_data,1) > size(mon_data,2)
        mon_data = mon_data';
    end
    n_of_models = size(mon_data, 1);
    output = zeros(n_of_models, ceil(size(mon_data, 2)/12));
    year_data_tmp = zeros(n_of_models, 1);
    year_ind = 1;
    month_count = 1;
    
    for i = 1: size(mon_data,2)

                year_data_tmp = year_data_tmp + mon_data(:, i);                 % summ of flow for a year
                if mod(month_count,12) == 0                                % when year is full:
                    output(:, year_ind) = year_data_tmp;                     % get result
                    year_data_tmp = zeros(n_of_models, 1);
                    year_ind = year_ind +1;
    %                 month_ind = 1;
                end 
                month_count = month_count +1;
    end
end

function[output] = filter_by_model(input, list_of_models)
    
    global list_of_models_hist
    output = input(logical(cat(2, 1, ismember(string(list_of_models_hist), string(list_of_models)))),:);
end

function[w_nv] = calc_w_nv(y_var, pc)
    if size(y_var) ~= size(pc)
        disp("size error! calc_w_nv function");
    end
    n_of_models = size(y_var, 1);
    a_tmp = zeros(n_of_models, 1);
    delta_a_tmp = a_tmp;
    for count = 1: n_of_models
        [a_tmp(count),delta_a_tmp(count)] = my_lin_reg(pc(count,:),y_var(count,:));
        w_nv(count) = exp(-((a_tmp(count)-a_tmp(1))^2)/(2*delta_a_tmp(1)^2));
    end  
end


function[output] = calc_y_f(y_var, pc)
    if size(y_var) ~= size(pc)
        disp("size error! calc_y_f function");
    end    
    n_of_models = size(y_var, 1);
    for count = 1: n_of_models
        [a_tmp,delta_a_tmp] = my_lin_reg(pc(count,:),y_var(count,:));
        output(count,:) = y_var(count,:) - mean(y_var(count,:)) - a_tmp*pc(count,:);
    end
end











































