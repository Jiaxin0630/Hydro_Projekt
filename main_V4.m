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
final.gravity_values = final.gravity_values - mean(final.gravity_values);

% ----------------------------------------------------------
figure(1),plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal with outliers')
legend("G1-F60","G2-F60")


%% Remove Outliers
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);
% ----------------------------------------------------------
figure(2), subplot(3,1,1),plot(final.time,diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors with outliers')
% ----------------------------------------------------------

[final.time,final.gravity_values] = rmoveoutliers(final.time,diffG1G2,final.gravity_values);

p = polyfit(seconds(final.time-final.time(1)),diffG1G2,1);
drift = seconds(final.time-final.time(1))*p(1);

% ----------------------------------------------------------
figure(2),subplot(3,1,2),plot(final.time, diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors without outliers')
% ----------------------------------------------------------

figure(2), subplot(3,1,3),plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal without outliers')
%% Remove drift
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);
final.gravity_values(:,2) = final.gravity_values(:,2) + drift;

%% Correcting tides 
idxx = length(final.gravity_values);

final.gravity_values(:,1)=final.gravity_values(:,1)- data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;
final.gravity_values(:,2)=final.gravity_values(:,2)-data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;

final.gravity_values = final.gravity_values - nanmean(final.gravity_values);

%% Correcting the influence of atmospheric pressure
data.Raw.Atmos.time = data.Raw.SG.time;
data.Raw.Atmos.cor = data.Raw.SG.PRESSURE_ADMIT_HPA_NMS2*1e-1*(data.Raw.SG.Br2_F60);
idx = data.Raw.Atmos.cor == 0;
final.gravity_values(idx,:) = NaN;
final.gravity_values = final.gravity_values - data.Raw.Atmos.cor;
final.gravity_values = final.gravity_values - nanmean(final.gravity_values);
legend("G1-F60","G2-F60")
% ----------------------------------------------------------
figure(5),
plot(final.time,final.gravity_values)
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
title('Gavity values in \muGal after correcting the influence of atmospheric pressure and tides (with remaining frequencies)')


%% Eliminate remaining frequencies
figure(6),subplot(2,1,1),plot(final.time,final.gravity_values(:,1)*10,'b--'),hold on,
subplot(2,1,2),plot(final.time,final.gravity_values(:,2)*10,'b--'),hold on,
f_remaing = [
%    1/(354*86400)
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
     1/(2*86400)
    %1/(13.1111*86400)
    ]; % day

G1_temp = final.gravity_values(~idx,1);
G2_temp = final.gravity_values(~idx,2);
[x1_hat,G1,Ax1] = FitSin(final.time(~idx),G1_temp,f_remaing);
[x2_hat,G2,Ax2] = FitSin(final.time(~idx),G2_temp,f_remaing);

subplot(2,1,1),plot(final.time(~idx),G1*10 )
title("G1-60")
legend("with remaining frequencies","without remaining frequencies")
subplot(2,1,2),plot(final.time(~idx),G2*10 )
title("G2-60")
legend("with remaining frequencies","without remaining frequencies")

final.time_temp = final.time(~idx);
final.gravity_values_temp = [G1 G2];

figure(7)
plot(final.time_temp,final.gravity_values_temp(:,1) - final.gravity_values_temp(1,1),'color',[0.8500 0.3250 0.0980],'LineWidth',2	)
hold on
plot(final.time_temp,final.gravity_values_temp(:,2) - final.gravity_values_temp(1,2),'color',[0.4660 0.6740 0.1880],'LineWidth',2)
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
title('Gavity values in \muGal after correcting the influence of atmospheric pressure and tides (without remaining frequencies)')
ylim([-4 8])
grid on
%% Soil Moisture
SMSolidus1 = loadSM("G:\Project_Hydro\BFO_data\Soil_Moisture\SMT100DataSolidus1.csv");
figure(10),
plot(data.Raw.Soil_Moisture.SMT100DataSolidus1.time ,SMSolidus1{1, 3}  )
ylabel('Soil Moisture [%]')
title("SMT100DataSolidus1")

SM = interp1(seconds(data.Raw.Soil_Moisture.SMT100DataSolidus1.time -...
    data.Raw.Soil_Moisture.SMT100DataSolidus1.time(1)),SMSolidus1{1, 3},seconds(final.time_temp-final.time_temp(1)));

%% Pre
len = length(data.Raw.Precipitation.pre)/60;
for i = 1:len-1
    pre_hour(i,1) = sum(data.Raw.Precipitation.pre((i-1)*60+1:(i)*60));
    time_hour(i,1) = data.Raw.Precipitation.time((i-1)*60+1);
end
figure
plot(time_hour,pre_hour,'LineWidth',2,'color','b')
ylabel('Precipitation in mm/hour')




%% Final
figure(8),
plot(final.time_alternative,final.gravity_alternative - mean(final.gravity_alternative))
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Gavity values in \muGal after correcting the influence of atmospheric pressure')
legend("G1-F60","G2-F60")


