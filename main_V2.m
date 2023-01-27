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
figure(1), subplot(2,1,1),plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal with outliers')
%% Remove Outliers
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);
figure, subplot(2,1,1),plot(final.time,diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors with outliers')

[final.time,final.gravity_values] = rmoveoutliers(final.time,diffG1G2,final.gravity_values);
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);


subplot(2,1,2),plot(final.time, diffG1G2)
ylabel('gravity values in [\muGal]')
title('Difference between two sensors without outliers')
figure(1), subplot(2,1,2),plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal without outliers')

%% Correcting tides
final.gravity_values(:,1)=final.gravity_values(:,1)-data.Raw.Synthetic_Tides.Tide(1:509760)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:509760)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:509760)/10;
final.gravity_values(:,2)=final.gravity_values(:,2)-data.Raw.Synthetic_Tides.Tide(1:509760)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:509760)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:509760)/10;

figure()
plot(final.time,final.gravity_values(:,1)-mean(final.gravity_values(:,1)))
hold on
plot(final.time,final.gravity_values(:,2)-mean(final.gravity_values(:,2)))
hold off
legend({'G1-F60','G2-F60'})
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
%axis([-inf inf -15 14])

%% The influence of atmospheric pressure

data.Raw.Atmos.time = data.Raw.SG.time;
data.Raw.Atmos.cor = data.Raw.SG.PRESSURE_ADMIT_HPA_NMS2*1e-1*(data.Raw.SG.Br2_F60);
data.Raw.Atmos.cor(data.Raw.Atmos.cor==0) = NaN;
ts = datetime(2022,1,1);
te = datetime(2022,12,1);
idx = data.Raw.Atmos.time>=ts & data.Raw.Atmos.time<te;
figure,set(gcf,"Units","normalized","Position",[0.26,0.28,0.48,0.3375]);
plot(data.Raw.Atmos.time(idx), data.Raw.Atmos.cor(idx),'LineWidth',2,'Color','b')
ylabel('Tide \muGal')
title('The influence of atmospheric pressure','FontSize',20)

%% correcting the influence of atmospheric pressure
data.Gravity_cal1Tides1Atmos1.time = data.Raw.SG.time;
final.gravity_values(:,1) = final.gravity_values(:,1) - data.Raw.Atmos.cor;
final.gravity_values(:,2) = final.gravity_values(:,2) - data.Raw.Atmos.cor;
idx=isnan(final.gravity_values(:,1));
%data.Gravity_cal1Tides1Atmos1.G1_F60 = data.Gravity_cal1Tides1.G1_F60 - data.Raw.Atmos.cor;
%data.Gravity_cal1Tides1Atmos1.G2_F60 = data.Gravity_cal1Tides1.G2_F60 - data.Raw.Atmos.cor;

%data_G1 = data.Gravity_cal1Tides1Atmos1.G1_F60;
%data_G2 = data.Gravity_cal1Tides1Atmos1.G2_F60;

figure()
plot(final.time(~idx),final.gravity_values(~idx,1)-mean(final.gravity_values(~idx,1)))
hold on
plot(final.time(~idx),final.gravity_values(~idx,2)-mean(final.gravity_values(~idx,2)))
hold off
legend({'G1-F60','G2-F60'})
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
axis tight
title('Calibrated gravity values in \muGal after removing atmospheric influence')
%axis([-inf inf -15 14])

%% Remaining tide frequencies
%periods
T = [1 0.5 1+50/60/24 (1+50/60/24)/2 167 334];
T=[1];
%figure()
% Y = fft(final.gravity_values(~idx,1));
% frequency = 1;
% T=1/frequency;
% L=length(final.gravity_values(~idx,1));
% t=(0:L-1)*T;
% P2=abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = frequency*(0:(L/2))/L;
% semilogy(f*60*24,P1)
% xlabel('frequency [cpd]')
% ylabel('power')
% axis([0 10 1e-3 1])
t = datenum(data.Gravity_cal1Tides1.time(~idx) )- datenum(data.Gravity_cal1Tides1.time(1) );
[P,f]=plomb(final.gravity_values(~idx,1),t);
figure()
semilogy(f,P)
axis([0 5 0 inf])
xlabel('frequency [cpd]')
ylabel('power spectral density')
for j =1:1
y=final.gravity_values(~idx,j);
for i = 1:length(T)
    phi = 0;
a = 0.3;
x=[phi; a];
M=zeros(length(y),1);
A=zeros(length(y),2);
while true
    for i = 1:length(y)
        M(i)=sin(2*pi*t(i)+phi);
        A(i,1)=cos(2*pi*t(i)+phi)*a;
        A(i,2)=sin(2*pi*t(i)+phi);
    end
    dy=y-M*a;
    dx=inv(A'*A)*A'*dy;
    xhat=x+dx;
    if(norm(dx)<1e-10)
        break
    else
        x;
        x=xhat;
        phi=x(1);
        a=x(2);
    end
end
yhat=M*a+A*dx;
ehat=y-yhat;
A'*ehat;
%     %linearized adjustment for the different frequencies
%     x=[0; 0.5];
%     M=sin(2*pi/T(i)*t+x(1));
%     A=[cos(2*pi/T(i)*t+x(1))*x(2) sin(2*pi/T(i)*t+x(1))];
%     while true
%         M=sin(2*pi/T(i)*t+x(1));
%         A=[cos(2*pi/T(i)*t+x(1))*x(2) sin(2*pi/T(i)*t+x(1))];
%         dy=y-M*x(2);
%         dx=(A'*A)\A'*dy;
%         xhat=x+dx;
%         if(norm(dx)<1e-8)
%             break
%         else
%             x=xhat;
%         end
%     end
%     yhat=M*x(2)+A*dx;
%     ehat=y-yhat;
%     y=ehat;
    %i
%     figure()
%     plot(final.time(~idx),y)
%     ylabel('gravity values in [\muGal]')
%     title('Calibrated gravity values in \muGal after eliminating tidal frequencies')
end
final.gravity_values(~idx,j)=y;
end
figure()
Y = fft(final.gravity_values(~idx,1));
frequency = 1;
T=1/frequency;
L=length(final.gravity_values(~idx,1));
t1=(0:L-1)*T;
P2=abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = frequency*(0:(L/2))/L;
semilogy(f*60*24,P1)
xlabel('frequency [cpd]')
ylabel('power')
axis([0 10 1e-3 1])
figure()
plot(final.time,final.gravity_values)
ylabel('gravity values in [\muGal]')
title('Calibrated gravity values in \muGal after eliminating tidal frequencies')
% A = [t...
%     cos(2*pi/T1*t) sin(2*pi/T1*t)...
%     cos(2*pi/T2*t) sin(2*pi/T2*t)...
%     cos(2*pi/T3*t) sin(2*pi/T3*t)...
%     cos(2*pi/T4*t) sin(2*pi/T4*t)];
% 
% par = (A'*A)\A'*data.Gravity_cal1Tides1.G1_F60;
%% drift
...




% figure,
% [data_G1,TFrm_G1,TFoutlier_G1] = rmoveoutliers(data_G1-nanmean(data_G1),"movmedian",hours(5),"SamplePoints",data.Gravity_cal1Tides1Atmos1.time);
% [data_G2,TFrm_G2,TFoutlier_G2] = rmoveoutliers(data_G2-nanmean(data_G2),"movmedian",hours(5),"SamplePoints",data.Gravity_cal1Tides1Atmos1.time);
% 
% plot(data.Gravity_cal1Tides1Atmos1.time(~TFrm_G1),data_G1)
% hold on,
% plot(data.Gravity_cal1Tides1Atmos1.time(~TFrm_G2),data_G2)
% 
% title('Gravity values over the whole year after subtraction of tides and atmospheric influence','FontSize',20)
% ylabel('gravity values [\muGal]')
% legend("G1-F60","G2-F60")
% 
% %% Analysis of soil moisture
% figure,set(gcf,'Units','normalized','Position',[0.14,0.11,0.74,0.72])
% fields = fieldnames(data.Raw.Soil_Moisture);
% for i = 1:numel(fields)
%     subplot(numel(fields),1,i)
%     station = fields{i};
%     IDX2 = data.Raw.Soil_Moisture.(station).time>=Ts & data.Raw.Soil_Moisture.(station).time<Te;
%     if ~contains(station,'TTGOPICO')
%         plot(data.Raw.Soil_Moisture.(station).time(IDX2),data.Raw.Soil_Moisture.(station).SM_Probe0(IDX2),'LineWidth',1.5,'Color','b')
%         hold on
%         plot(data.Raw.Soil_Moisture.(station).time(IDX2),data.Raw.Soil_Moisture.(station).SM_Probe1(IDX2),'LineWidth',1.5,'Color','r')
%         legend('SM Probe0','SM Probe1')
%     else
%         plot(data.Raw.Soil_Moisture.(station).time(IDX2),data.Raw.Soil_Moisture.(station).SM(IDX2),'LineWidth',1.5,'Color','b')
%         legend('SM')
%     end
%     try
%         
%     catch
% 
%     end
%     ylabel('SM in [%]')
%     title(['Station: ' station],'FontSize',20)
% end
% %% Analysis of prepcipitation
% figure,set(gcf,'Units','normalized','Position',[0.23,0.44,0.55,0.32])
% plot(data.Raw.Precipitation.time,data.Raw.Precipitation.pre,'LineWidth',2,'Color','b')
% 
% ylabel('precipitation [mm/min]')
% title('The precipitation gauge','FontSize',20)
% %% plot
% 
% figure,set(gcf,'Units','normalized','Position',[0,0,1,0.93])
% subplot(6,1,1),plot(data.Raw.SG.time,data.Raw.SG.G1_F60,'LineWidth',2,'Color','b')
% hold on,plot(data.Gravity_cal1.time,data.Raw.SG.G2_F60,'LineWidth',2,'Color','r')
% 
% legend('G1-F60','G2-F60')
% ylabel('gravity values in volts from sensor [v]')
% title('Uncalibrated gravity values','FontSize',20)
% 
% subplot(6,1,2),plot(data.Gravity_cal1.time,data.Gravity_cal1.G1_F60,'LineWidth',2,'Color','b')
% hold on,plot(data.Gravity_cal1.time,data.Gravity_cal1.G2_F60,'LineWidth',2,'Color','r')
% 
% legend('G1-F60','G2-F60')
% ylabel('gravity values [\muGal]')
% title('Calibrated gravity values','FontSize',20)
% 
% 
% subplot(6,1,3),plot(data.Gravity_cal1Tides1.time,data.Gravity_cal1Tides1.G1_F60-mean(data.Gravity_cal1Tides1.G1_F60),'LineWidth',2,'Color','b')
% hold on,plot(data.Gravity_cal1Tides1.time,data.Gravity_cal1Tides1.G2_F60-mean(data.Gravity_cal1Tides1.G2_F60),'LineWidth',2,'Color','r')
% 
% legend('G1-F60','G2-F60')
% ylabel('gravity values [\muGal]')
% title('Gravity values with the effect of tides removed','FontSize',20)
% 
% 
% subplot(6,1,4),plot(data.Gravity_cal1Tides1Atmos1.time,data.Gravity_cal1Tides1Atmos1.G1_F60 -nanmean(data.Gravity_cal1Tides1Atmos1.G1_F60 ),'LineWidth',2,'Color','b')
% hold on,plot(data.Gravity_cal1Tides1Atmos1.time,data.Gravity_cal1Tides1Atmos1.G2_F60 -nanmean(data.Gravity_cal1Tides1Atmos1.G2_F60 ),'LineWidth',2,'Color','r')
% 
% legend('G1-F60','G2-F60')
% ylabel('gravity values [\muGal]')
% title('Gravity values with the effect of tides, atmospheric pressure and mean removed','FontSize',20)
% 
% station = 'SMT100DataSolidus1'; % This station is located near SG
% IDX2 = data.Raw.Soil_Moisture.(station).time>=Ts & data.Raw.Soil_Moisture.(station).time<Te;
% subplot(6,1,5),plot(data.Raw.Soil_Moisture.(station).time(IDX2),...
%     data.Raw.Soil_Moisture.(station).SM_Probe0(IDX2),'LineWidth',2,'Color','b')
% hold on,plot(data.Raw.Soil_Moisture.(station).time(IDX2),...
%     data.Raw.Soil_Moisture.(station).SM_Probe1(IDX2),'LineWidth',2,'Color','r')
% try
%     
%     legend("SM Probe0","SM Probe1")
% catch
%     
% end
% ylabel('Soil moisture in [%]')
% title('Soil moisture in station ImkoDataSolidus5','FontSize',20)
% 
% subplot(6,1,6),plot(data.Raw.Precipitation.time,data.Raw.Precipitation.pre,'LineWidth',2,'Color','b')
% 
% ylabel('precipitation [mm/min]')
% title('The precipitation gauge','FontSize',20)
% 
% 
% figure,plot(data.Gravity_cal1Tides1Atmos1.time,data.Gravity_cal1Tides1Atmos1.G1_F60 -mean(data.Gravity_cal1Tides1Atmos1.G1_F60 ),'LineWidth',2,'Color','b')
% %% Plot
% %% Uncalibrated and calibrated Gravity values from sensor G1-F60 and G2-F60
% figure(1),set(gcf,'unit','normalized','Position',[0.25,0.23,0.47,0.5])
% subplot(2,1,1),plot(data.Gravity_cal1.time,data.Raw.SG.G1_F60,'LineWidth',2,'Color','b')
% hold on,plot(data.Raw.SG.time,data.Raw.SG.G2_F60,'LineWidth',2,'Color','r')
% legend('G1-F60: uncalibrated','G2-F60: uncalibrated')
% ylabel('gravity values in volts from sensor [v]')
% title('Uncalibrated gravity values in volts from sensor','FontSize',20)
% 
% subplot(2,1,2),plot(data.Gravity_cal1.time,data.Gravity_cal1.G1_F60,'LineWidth',2,'Color','b')
% hold on,plot(data.Gravity_cal1.time,data.Gravity_cal1.G2_F60,'LineWidth',2,'Color','r')
% 
% legend('G1-F60: Calibrated','G2-F60: Calibrated')
% ylabel('gravity values in volts from sensor [\muGal]')
% title('Calibrated gravity values in \muGal from sensor','FontSize',20)
% %% 
