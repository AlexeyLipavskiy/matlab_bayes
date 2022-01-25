clear all
close all
clc
%%

load hist+ssp119+245+585.mat

global w_length
w_length = 1000; % length of distribution in points
%% just a plot
figure;
h = plot(years,P_245,'LineWidth',2);
legend(list_of_models_245)

% h(1).LineWidth = 2.5;
% все веса считать до 14 года, то есть по данным хисторикал. дельта ка
% (каэф лин тренда и его стандартное отклонение делать по книге
% НумРекФортран стр 686(там другие обозначения). Чуть позже построить корреляционные функции
% сигнала и аппроксимировать экспонентой, найти значение "декремента" -
% когда амплитуда падает в е раз, потом пересчитать дельта к по формуле в
% книге (там будет суммирование не по всему периоду, а по 2*декремент.

%% removing Nan rows

% P_119(any(isnan(P_119(:,1)), 2), :) = [];
% R_119(any(isnan(R_119(:,1)), 2), :) = [];
% Rs_119(any(isnan(Rs_119(:,1)), 2), :) = [];
% 
% P_245(any(isnan(P_245(:,1)), 2), :) = [];
% R_245(any(isnan(R_245(:,1)), 2), :) = [];
% Rs_245(any(isnan(Rs_245(:,1)), 2), :) = [];
% 
% P_585(any(isnan(P_585(:,1)), 2), :) = [];
% R_585(any(isnan(R_585(:,1)), 2), :) = [];
% Rs_585(any(isnan(Rs_585(:,1)), 2), :) = [];
%% fix years
years = years(1:37);

P_119 = P_119(:,1:37);
R_119 = R_119(:,1:37);
Rs_119 = Rs_119(:,1:37);

P_245 = P_245(:,1:37);
R_245 = R_245(:,1:37);
Rs_245 = Rs_245(:,1:37);

P_585 = P_585(:,1:37);
R_585 = R_585(:,1:37);
Rs_585 = Rs_585(:,1:37);


%% 

% reg = fitlm(years, P_585(1,:))
%%



[P_585_w_m,P_585_w_tr,P_585_w_iav] = calc_all_models(years, P_585);





%%

plot(P_585_w_m(1,:));

%%

% figure
% plot(normpdf(1:w_length,P_585_m,P_585_std))
% a = normpdf(1:w_length,P_585_m,P_585_std);

%%

% tmp = P_585(1,:);
% cor_tmp = ifft(fft(tmp).*conj(fft(tmp)));
% PlotCorr(cor,1);

% 
% figure;
% plot_data_and_trend(P_585(1,:));

% count = 1;
% data = P_585;
% yrs = years;




function [data_w_m,data_w_tr,data_w_iav] = calc_all_models(yrs, data)
num_of_models = size(data,1);
global w_length
data_mean = zeros(num_of_models,w_length);
data_std = zeros(num_of_models,w_length);
data_k = zeros(num_of_models,w_length);
data_dk = zeros(num_of_models,w_length);
data_dy = zeros(num_of_models,w_length);
data_ystd = zeros(num_of_models,w_length);
data_w_m = zeros(num_of_models,w_length);
data_w_tr = zeros(num_of_models,w_length);
data_w_iav = zeros(num_of_models,w_length);

    for count = 1 : num_of_models

%         [data_m(model_num),data_std(model_num),data_k(model_num),data_dk(model_num),data_dy(model_num),data_ystd(model_num),...
%             data_w_m(model_num),data_w_tr(model_num),data_w_iav(model_num)] = my_lin_reg(yrs, data(model_num,:));
        [data_w_m_tmp,data_w_tr_tmp,data_w_iav_tmp] = my_lin_reg(yrs, data(count,:));

        data_w_m(count,:) = data_w_m_tmp;
        data_w_tr(count,:) = data_w_tr_tmp;
        data_w_iav(count,:) = data_w_iav_tmp;
    end
end



function [w_m,w_tr,w_iav] = my_lin_reg(x,y)
%           Y=A+B*X
% Using symbols as in the book "NumRecFortran" page 655
global w_length
S_x = mean(x,2);
S_y = mean(y,2);

ndata = length(x);

sx = 0;
sy = 0;
ss = 1;
st2 = 0;
b = 0;
for count = 1 : ndata
    sx = sx +x(count);
    sy = sy+ y(count);
end
ss = ndata;
sxoss = sx/ss;
for count = 1 : ndata
    t = x(count) - sxoss;
    st2 = st2 + t^2;
    b = b + t*y(count);  
end
b = b/st2;
a = (sy-sx*b)/ss;
siga = sqrt((1 + sx*sx/(ss*st2))/ss);
sigb = sqrt(1/st2);

y_mean = S_y;  
y_std = std(y,0,2)
k = b;
delta_k = sigb;
delta_y = detrend(y);
delta_y_std = std(delta_y,0,2);

w_m = normpdf(y_mean-3*y_std:6*y_std/w_length:y_mean+3*y_std-6*y_std/w_length,y_mean,y_std);
w_tr = normpdf(k-3*delta_k:6*delta_k/w_length:k+3*delta_k-6*delta_k/w_length,k,delta_k);
w_iav = NaN;
end


function plot_data_and_trend (data)
plot([1:length(data)],polyval(polyfit([1:length(data)]',data,1),[1:length(data)]),[1:length(data)],data);
end


function [y_mean,y_std,regcoef,regcoefstd,regr,regss] = lin_reg(x,y)
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
  coefdof = (1-abs(rx*ry))/(1+abs(rx*ry));
  eqn = length(x)*coefdof;
  regcoefstd = sqrt((y_variance/x_variance-regcoef^2)/eqn);
  dof = (length(x)-2)*coefdof;
  regss = 0;
  if(isnan(regcoefstd))
    [regcoefstd x_mean y_mean x_variance y_variance]
  end
end
































