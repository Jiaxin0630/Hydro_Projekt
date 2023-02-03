clear all
close all
clc
%addpath BFO_data
load data2.mat

%% Initialize the final result
% ----------------------------------------------------------
% ----------------------------------------------------------
% ----------------------------------------------------------
final.gravity_values = [];
final.time = [];

% ----------------------------------------------------------
% ----------------------------------------------------------
% ----------------------------------------------------------
%% Uncalibrated and calibrated Gravity values from sensor G1-F60 and G2-F60

data.Gravity_cal1.time = data.Raw.SG.time;
data.Gravity_cal1.G1_F60 = data.Raw.SG.G1_F60*data.Raw.SG.GRAVITY_CAL_1_UGAL_V;
data.Gravity_cal1.G2_F60 = data.Raw.SG.G2_F60*data.Raw.SG.GRAVITY_CAL_2_UGAL_V;

final.time = data.Gravity_cal1.time;
final.gravity_values = [data.Gravity_cal1.G1_F60 data.Gravity_cal1.G2_F60];
% final.gravity_values = final.gravity_values - final.gravity_values(1,:);
final.gravity_values = final.gravity_values - mean(final.gravity_values);

% ----------------------------------------------------------
figure(1), subplot(6,1,1),
plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal with outliers')
legend('G1-F60','G2-F60')
figure(7),subplot(6,1,1),
plot(final.time,final.gravity_values - mean(final.gravity_values))
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal with outliers')
legend('G1-F60','G2-F60')
% ----------------------------------------------------------

%% Remove Outliers
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);
% ----------------------------------------------------------
figure(2), subplot(2,1,1),
plot(final.time,diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors with outliers')
% ----------------------------------------------------------

[final.time,final.gravity_values] = rmoveoutliers(final.time,diffG1G2,final.gravity_values);
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);

% ----------------------------------------------------------
figure(2), subplot(2,1,2),
plot(final.time, diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors without outliers')
% ----------------------------------------------------------

figure(1), subplot(6,1,2),plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal without outliers')
legend('G1-F60','G2-F60')
figure(7),subplot(6,1,2),
plot(final.time,final.gravity_values - mean(final.gravity_values))
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal without outliers')
legend('G1-F60','G2-F60')
% ----------------------------------------------------------
%% 



%% Correcting tides 
idxx = length(final.gravity_values);

final.gravity_values(:,1)=final.gravity_values(:,1)-data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;
final.gravity_values(:,2)=final.gravity_values(:,2)-data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;

% ---------------------------------------------------------
data.Raw.Atmos.time = data.Raw.SG.time;
data.Raw.Atmos.cor = data.Raw.SG.PRESSURE_ADMIT_HPA_NMS2*1e-1*(data.Raw.SG.Br2_F60);
idx = data.Raw.Atmos.cor==0;
final.gravity_values(idx) = NaN;
final.gravity_values = final.gravity_values - data.Raw.Atmos.cor;

% ----------------------------------------------------------
figure(1),subplot(6,1,3)
plot(final.time,final.gravity_values(:,1))
hold on
plot(final.time,final.gravity_values(:,2))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting tides with remaining frequencies')

figure(7),subplot(6,1,3)
plot(final.time,final.gravity_values(:,1)-mean(final.gravity_values(:,1)))
hold on
plot(final.time,final.gravity_values(:,2)-mean(final.gravity_values(:,2)))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting tides with remaining frequencies')

%axis([-inf inf -15 14])
% ----------------------------------------------------------
t=datenum(final.time(~idx));
plomb(final.gravity_values(~idx,1),t)
%N = length(final.gravity_values(:,1));
% P2 = abs(H/N);P1 = P2(1:N/2+1);P1(2:end-1) = 2*P1(2:end-1);
% f = 1/60*(0:(N/2))/N;
% figure(3),
% % subplot(2,1,1),semilogy(1./f/86400,P1,'-o','LineWidth',1);
% subplot(2,1,1),semilogy(1./f/86400,P1,'-o','LineWidth',1);
% xlabel('Period [day]'), ylabel('Amplitude')
% title('The main period of the the gavity values in \muGal after correcting tides with remaining frequencies')
% legend('G1-F60')

% [P1,f] = periodogram(final.gravity_values(:,1),[],[],1/60);

% [~,idx] = maxk(P1,5);
% f_remaing = f(idx(2:end));
% day_remaing = 1./(f(idx(2:end)))/86400;

%% Eliminate remaining frequencies
figure(4),subplot(2,1,1),plot(final.time,final.gravity_values(:,1),'b--')
subplot(2,1,2),plot(final.time,final.gravity_values(:,2),'b--')
f_remaing = [
    %1/(354*86400)
%     1/(177*86400)
%     1/(118*86400)
%     1/(88.5003*86400)
%     1/(59*86400)
%     1/(19.6667*86400)
%     1/(3.92341*86400)
     1/(1*86400+50*60)
    1/86400
    1/(0.5*86400)
    2/(1*86400+50*60)
    % 1/(2*86400)

    %1/(13.1111*86400)
    ]; % day

final.gravity_values = final.gravity_values - nanmean(final.gravity_values);
G1 = final.gravity_values(:,1);
G2 = final.gravity_values(:,2);
figure,plot(final.time(~idx),G1(~idx))
[x1_hat,G11,Ax1] = FitSin(final.time(~idx),G1(~idx),f_remaing);
[x2_hat,G22,Ax2] = FitSin(final.time(~idx),G2(~idx),f_remaing);

hold on,plot(final.time(~idx),Ax1)

figure
t=datenum(final.time(~idx));
plomb(G11,t);


figure(4),subplot(2,1,1),hold on,plot(final.time,final.gravity_values(:,1),'r','LineWidth',1.5)
title('G1-F60')
legend('with remaining frequencies','without remaining frequencies')
subplot(2,1,2),hold on,plot(final.time,final.gravity_values(:,2),'r','LineWidth',1.5)
title('G2-F60')
legend('with remaining frequencies','without remaining frequencies')

H = fft(final.gravity_values(:,1));N = length(final.gravity_values(:,1));
P2 = abs(H/N);P1 = P2(1:N/2+1);P1(2:end-1) = 2*P1(2:end-1);
f = 1/60*(0:(N/2))/N;
figure(3),
%subplot(2,1,2),semilogy(1./f/86400,P1,'-o','LineWidth',1);
subplot(2,1,2),semilogy(f,P1,'-o','LineWidth',1);
xlabel('Period [day]'), ylabel('Amplitude')
title('The main period of the the gavity values in \muGal after correcting tides without remaining frequencies')
legend('G1-F60')

figure(1),subplot(6,1,4)
plot(final.time,final.gravity_values(:,1))
hold on
plot(final.time,final.gravity_values(:,2))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting tides without remaining frequencies')

figure(7),subplot(6,1,4)
plot(final.time,final.gravity_values(:,1)-mean(final.gravity_values(:,1)))
hold on
plot(final.time,final.gravity_values(:,2)-mean(final.gravity_values(:,2)))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting tides without remaining frequencies')

%% remove trend
figure(5),
subplot(2,1,1),plot(final.time,final.gravity_values(:,1),'b--'),hold on
subplot(2,1,2),plot(final.time,final.gravity_values(:,2),'b--'),hold on

d1 = detrend(final.gravity_values(:,1),1);
d2 = detrend(final.gravity_values(:,2),1);

final.gravity_values(:,1) = d1 - (d1(1) - final.gravity_values(1,1));
final.gravity_values(:,2) = d2 - (d2(1) - final.gravity_values(1,2));

subplot(2,1,1),plot(final.time,final.gravity_values(:,1),'r','LineWidth',2)
title('G1-F60')
subplot(2,1,2),plot(final.time,final.gravity_values(:,2),'r','LineWidth',2)
title('G2-F60')

figure(1),subplot(6,1,5)
plot(final.time,final.gravity_values(:,1))
hold on
plot(final.time,final.gravity_values(:,2))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after detrending')
figure(7),subplot(6,1,5)
plot(final.time,final.gravity_values - mean(final.gravity_values));
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after detrending')
%% Find own atmospheric factor
atmos = data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60 ~= 0);
final.gravity_alternative = final.gravity_values(data.Raw.SG.Br2_F60 ~= 0,:);
final.time_alternative = final.time(data.Raw.SG.Br2_F60 ~= 0,:);

figure(6),subplot(1,2,1)
x = atmos - mean(atmos);
y1 = final.gravity_alternative(:,1) - mean(final.gravity_alternative(:,1));
density1 = dscatter(x,y1,'MARKER',"o"); colorbar
xlabel('pressure - mean [hPa]')
ylabel('gravity - mean [\muGal]')
title('G1-F60')


figure(6),subplot(1,2,2)
y2 = final.gravity_alternative(:,2) - mean(final.gravity_alternative(:,2));
density2 = dscatter(x,y2,'MARKER',"o"); colorbar
xlabel('pressure - mean [hPa]')
ylabel('gravity - mean [\muGal]')
title('G2-F60')

P1 = speye(length(density1));
P1 = spdiags(1./density1,0,P1);

atmos_fac1 = FitLine(x,y1,P1);

P2 = speye(length(density2));
P2 = spdiags(1./density2,0,P2);

atmos_fac2 = FitLine(x,y2,P2);

%% correcting the influence of atmospheric pressure
atmos_G1 = atmos_fac1*atmos;
atmos_G2 = atmos_fac2*atmos;

final.gravity_alternative = final.gravity_alternative + [atmos_G1 atmos_G2];
figure(1), subplot(6,1,6),
plot(final.time_alternative,final.gravity_alternative)
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting the influence of atmospheric pressure')
figure(7), subplot(6,1,6),
plot(final.time_alternative,final.gravity_alternative - mean(final.gravity_alternative))
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting the influence of atmospheric pressure')
%% Final
figure(8),
plot(final.time_alternative,final.gravity_alternative - mean(final.gravity_alternative))
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting the influence of atmospheric pressure')



