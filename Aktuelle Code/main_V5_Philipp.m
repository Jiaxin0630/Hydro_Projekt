clear all
close all
clc
%addpath C:\Users\Jiaxin\Desktop\Hydro_Projekt
addpath SHbundle-master
%addpath uberall-master
%addpath BFO_data/Soil_Moisture
load data2.mat

plotting = 0;
later = 1;

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
%final.gravity_values = final.gravity_values - mean(final.gravity_values);

%plotting uncalibrated values
if plotting
figure(20)
subplot(2,1,1)
plot(final.time,data.Raw.SG.G1_F60)
ylabel('gravity meas. [V]')
legend("G1-F60")
subplot(2,1,2)
plot(final.time,data.Raw.SG.G2_F60)
ylabel('gravity meas. [V]')
legend("G2-F60")
end

% ----------------------------------------------------------
if plotting
%plot of calibrated gravity series fo full year
    figure(1),plot(final.time,final.gravity_values)
ylabel('gravity values [\muGal]')
%title('Calibrated gravity values in \muGal with outliers')
legend("G1-F60","G2-F60")

%plot for April 8th
figure(21),plot(final.time(139681:141120),final.gravity_values(139681:141120,:))
ylabel('gravity values [\muGal]')
%title('Calibrated gravity values in \muGal with outliers')
legend("G1-F60","G2-F60")
end

% if plotting
% figure(20)
% subplot(2,1,1)
% plot(final.time,final.gravity_values(:,1))
% ylabel('gravity meas. [\muGal]')
% legend("G1-F60")
% subplot(2,1,2)
% plot(final.time,final.gravity_values(:,2))
% ylabel('gravity meas. [\muGal]')
% legend("G2-F60")
% end

%% Remove Outliers
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);

if plotting
figure(22)
plot(final.time,diffG1G2)
ylabel('gravity difference [\muGal]')
%title('Differences between two sensors')
% ----------------------------------------------------------
figure(2), subplot(3,1,1),plot(final.time,diffG1G2)
ylabel('grav. difference [\muGal]')
title('Difference between two sensors with outliers')
% ----------------------------------------------------------
figure(23)
subplot(2,1,1)
plot(final.time,final.gravity_values)
ylabel('grav. values [\muGal]')
title('Gravity series before outlier removal')
legend(["G1-F60","G2-F60"],'Location','southeast')
end
[final.time,final.gravity_values] = rmoveoutliers(final.time,diffG1G2,final.gravity_values);

% ----------------------------------------------------------
if plotting
diffG1G2 = final.gravity_values(:,1) - final.gravity_values(:,2);
figure(2),subplot(3,1,2),plot(final.time, diffG1G2)
ylabel('grav. diff. [\muGal]')
title('Difference between two sensors without outliers')
% ----------------------------------------------------------

figure(2), subplot(3,1,3),plot(final.time,final.gravity_values)
ylabel('grav. values [\muGal]')
title('Calibrated gravity values in \muGal without outliers')

figure(23)
subplot(2,1,2)
plot(final.time,final.gravity_values)
ylabel('grav. values [\muGal]')
title('Gravity series after outlier removal')
legend(["G1-F60","G2-F60"],'Location','southeast')
end

%% Correcting tides 
idxx = length(final.gravity_values);

if plotting
%comparison before and after tide correction, full year
figure(25)
subplot(2,1,1)
plot(final.time,final.gravity_values)
ylabel('grav. values [\muGal]')
title('Gravity series before tide subtraction')
legend(["G1-F60","G2-F60"],'Location','southeast')

%comparison before and after tide correction, April 8th
figure(26)
subplot(2,1,1)
yyaxis left
plot(final.time(139681:141120),final.gravity_values(139681:141120,1))
ylabel('grav. values [\muGal]')
ylim([50 200])
yyaxis right
plot(final.time(139681:141120),final.gravity_values(139681:141120,2))
ylabel('grav. values [\muGal]')
ylim([10 160])
title('Gravity series before tide subtraction')
legend(["G1-F60","G2-F60"],'Location','southeast')
end

final.gravity_values(:,1)=final.gravity_values(:,1)- data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;
final.gravity_values(:,2)=final.gravity_values(:,2)-data.Raw.Synthetic_Tides.Tide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10-...
    data.Raw.Synthetic_Tides.LODTide(1:idxx)/10;

%plotting tides
if plotting
figure(24)
subplot(3,1,1)
plot(final.time,data.Raw.Synthetic_Tides.Tide(1:idxx)/10)
ylabel('tides [\muGal]')
%title('Gravity series after outlier removal')
%legend(["G1-F60","G2-F60"],'Location','southeast')
subplot(3,1,2)
plot(final.time,data.Raw.Synthetic_Tides.PoleTide(1:idxx)/10)
ylabel('Pole tides [\muGal]')
subplot(3,1,3)
plot(final.time,data.Raw.Synthetic_Tides.LODTide(1:idxx)/10)
ylabel('LOD tides [\muGal]')

figure(25)
subplot(2,1,2)
yyaxis left
plot(final.time,final.gravity_values(:,1))
ylim([102.5 127.5])
ylabel('grav. values [\muGal]')
title('Gravity series after tide subtraction')
yyaxis right
plot(final.time,final.gravity_values(:,2))
ylim([67.5 92.5])
ylabel('grav. values [\muGal]')
legend(["G1-F60","G2-F60"],'Location','southeast')

%comparison before and after tide correction, April 8th
figure(26)
subplot(2,1,2)
yyaxis left
plot(final.time(139681:141120),final.gravity_values(139681:141120,1))
ylabel('grav. values [\muGal]')
ylim([117 121])
yyaxis right
plot(final.time(139681:141120),final.gravity_values(139681:141120,2))
ylabel('grav. values [\muGal]')
ylim([81 85])
title('Gravity series after tide subtraction')
legend(["G1-F60","G2-F60"],'Location','southeast')
end

%final.gravity_values = final.gravity_values - nanmean(final.gravity_values);

%% AR
mp=mean(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900));
if plotting
figure(16)
scatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,final.gravity_values(data.Raw.SG.Br2_F60>900,1),1,'filled')
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean')
%axis([-inf inf -7 7])
figure(17)
scatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,final.gravity_values(data.Raw.SG.Br2_F60>900,2),1,'filled')
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean')
%axis([-inf inf -7 7])
end
idx = data.Raw.SG.Br2_F60 == 0;
t=datenum(final.time)-datenum(final.time(1));
if plotting
%adjustments without drift removed
figure(18)
co1=dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,final.gravity_values(data.Raw.SG.Br2_F60>900,1)-mean(final.gravity_values(data.Raw.SG.Br2_F60>900,1)));
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean, sensor 1')
ylim([-11 11])

figure(19)
co2 = dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,final.gravity_values(data.Raw.SG.Br2_F60>900,2)-mean(final.gravity_values(data.Raw.SG.Br2_F60>900,2)));
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean, sensor 2')
ylim([-11 11])
%axis([-inf inf -7 7])

Co1 = speye(length(co1));Co2 = speye(length(co2));
P1 = spdiags(co1,0,Co1);P2 = spdiags(co2,0,Co2);

A=data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-mp;
y1=final.gravity_values(data.Raw.SG.Br2_F60>900,1)-mean(final.gravity_values(data.Raw.SG.Br2_F60>900,1));
y2=final.gravity_values(data.Raw.SG.Br2_F60>900,2)-mean(final.gravity_values(data.Raw.SG.Br2_F60>900,2));
par1 = inv(A'*P1*A)*A'*P1*y1
par2 = inv(A'*P1*A)*A'*P2*y2

figure(18)
text(-24,-9.5,strcat('atmospheric factor:', {' '}, num2str(par1), ' \muGal/hPa'),'FontSize',10)
figure(19)
text(-24,-9.5,strcat('atmospheric factor:', {' '}, num2str(par2), ' \muGal/hPa'),'FontSize',10)

%adjustments with drift removed
p1=polyfit(t(~idx),final.gravity_values(~idx,1),1);
grav1=final.gravity_values(~idx,1)-polyval(p1,t(~idx));
p2=polyfit(t(~idx),final.gravity_values(~idx,2),1);
grav2=final.gravity_values(~idx,2)-polyval(p2,t(~idx));

figure(35)
co1=dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,grav1);
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean, sensor 1')
ylim([-8 8])

figure(36)
co2 = dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-...
    mp,grav2);
xlabel('pressure [hPa] - mean')
ylabel('gravity [\muGal] - mean, sensor 2')
ylim([-8 8])
%axis([-inf inf -7 7])

Co1 = speye(length(co1));Co2 = speye(length(co2));
P1 = spdiags(co1,0,Co1);P2 = spdiags(co2,0,Co2);

A=data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)-mp;
y1=grav1;
y2=grav2;
par1 = inv(A'*P1*A)*A'*P1*y1
par2 = inv(A'*P1*A)*A'*P2*y2
figure(35)
text(-24,-6,strcat('atmospheric factor:', {' '}, num2str(par1), ' \muGal/hPa'),'FontSize',10)
figure(36)
text(-24,-6,strcat('atmospheric factor:', {' '}, num2str(par2), ' \muGal/hPa'),'FontSize',10)
end
%% Correcting the influence of atmospheric pressure
if plotting
%plot pressure series
figure(27)
plot(final.time,data.Raw.SG.Br2_F60)
ylabel('pressure [hPa]')

%finding alternative for outside pressure: inside pressure
figure(28)
dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900),data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900));
xlabel('outside pressure [hPa]')
ylabel('inside pressure [hPa]')
c=corrcoef(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900),data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900));
text(927,968,strcat('Corr. coefficient:', {' '}, num2str(c(1,2))),'FontSize',9)

%finding alternative for outside pressure: pressure at magnetic hut
data.Raw.MagneticHutPressure = loadMagneticHut(final.time);
figure(29)
plot(final.time,data.Raw.MagneticHutPressure)
ylabel('pressure [hPa]')
end

%comparison outside/inside pressure and gravity residuals
figure(30)
subplot(2,1,1)
plot(final.time,[data.Raw.SG.Br1_F60 data.Raw.SG.Br2_F60])
legend('inside pressure','outside pressure')
ylabel('pressure [hPa]')
ylim([920 980])
subplot(2,1,2)
plot(final.time,final.gravity_values)
legend('G1-F60','G2-F60')
ylabel('gravity residual [\muGal]')
%ylim([920 980])

%which pressure series match gravity residuals better?
figure(31)
dscatter(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,1));
xlabel('outside pressure [hPa]')
ylabel('gravity resiuals [\muGal], sensor 1')
c=corrcoef(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,1));
text(927,108,strcat('Corr. coefficient:', {' '}, num2str(c(1,2))),'FontSize',9)
%corrcoef(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,2))

%which pressure series match gravity residuals better?
figure(32)
dscatter(data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,1));
xlabel('inside pressure [hPa]')
ylabel('gravity resiuals [\muGal], sensor 1')
c=corrcoef(data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,1));
text(927,108,strcat('Corr. coefficient:', {' '}, num2str(c(1,2))),'FontSize',9)
%corrcoef(data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900),final.gravity_values(data.Raw.SG.Br2_F60>900,2))

%plot atmospheric correction on April 8th
figure(34)
subplot(3,1,1)
yyaxis left
plot(final.time(139681:141120),final.gravity_values(139681:141120,1))
ylabel('gravity res. [\muGal]')
ylim([117 121])
yyaxis right
plot(final.time(139681:141120),final.gravity_values(139681:141120,2))
ylabel('gravity res. [\muGal]')
ylim([81 85])
legend({"G1-F60","G2-F60"},'Location','southeast')


%correction
data.Raw.Atmos.time = data.Raw.SG.time;
data.Raw.Atmos.cor = data.Raw.SG.PRESSURE_ADMIT_HPA_NMS2*1e-1*(data.Raw.SG.Br2_F60-mean(data.Raw.SG.Br2_F60(data.Raw.SG.Br2_F60>900)));
data.Raw.Atmos.cor1 = data.Raw.Atmos.cor;
data.Raw.Atmos.cor1(idx) = data.Raw.SG.PRESSURE_ADMIT_HPA_NMS2*1e-1*(data.Raw.SG.Br1_F60(idx)-mean(data.Raw.SG.Br1_F60(data.Raw.SG.Br2_F60>900)));
final.gravity_values1 = final.gravity_values - data.Raw.Atmos.cor1;
final.gravity_values(idx,:) = NaN;
final.gravity_values = final.gravity_values - data.Raw.Atmos.cor;
final.gravity_values = final.gravity_values - nanmean(final.gravity_values);
final.gravity_values1 = final.gravity_values1 - mean(final.gravity_values1);

%plot with inside pressure corrected
if plotting
figure(33)
plot(final.time,final.gravity_values1(:,1)-mean(final.gravity_values1(:,1)))
hold on
plot(final.time,final.gravity_values1(:,2)-mean(final.gravity_values1(:,2)))
hold off
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')

%plot April 8th
figure(34)
subplot(3,1,2)
plot(final.time(139681:141120),data.Raw.SG.Br2_F60(139681:141120))
ylabel('outside press. [hPa]')

subplot(3,1,3)
yyaxis left
plot(final.time(139681:141120),final.gravity_values1(139681:141120,1))
ylabel('gravity res. [\muGal]')
ylim([-2.75 -2.1])
yyaxis right
plot(final.time(139681:141120),final.gravity_values1(139681:141120,2))
ylabel('gravity res. [\muGal]')
ylim([-5.18 -4.53])
legend('G1-F60','G2-F60')
end
% ----------------------------------------------------------
if plotting
figure(5),
plot(final.time,final.gravity_values)
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
title('Gravity values in \muGal after correcting the influence of atmospheric pressure and tides (with remaining frequencies)')
end


%% Eliminate remaining frequencies
if plotting
%search for frequencies
[P,f]=plomb(final.gravity_values1(:,1),t);
figure(37)
semilogy(f,P)
axis([0 5 1e-5 inf])
%title('')
xlabel('frequency [cpd]')
ylabel('power spectral density')

figure(44)
subplot(2,1,1)
semilogy(f,P)
axis([0 5 1e-7 inf])
title('PSD before eliminating frequencies')
xlabel('frequency [cpd]')
ylabel('power spectral density')

figure(39)
loglog(1./f,P)
xlabel('period [d]')
ylabel('power spectral density')
axis tight

%sensor 2
[P,f]=plomb(final.gravity_values1(:,2),t);
figure(40)
semilogy(f,P)
axis([0 5 1e-5 inf])
%title('')
xlabel('frequency [cpd]')
ylabel('power spectral density')

figure(41)
loglog(1./f,P)
xlabel('period [d]')
ylabel('power spectral density')
axis tight

figure(38)
plot(final.time(129601:129600+30*60*24),final.gravity_values1(129601:129600+30*60*24,:))
ylabel('gravity residuals [\muGal]')
legend('G1-F60','G2-F60')

%eliminate frequencies
figure(6),subplot(2,1,1),plot(final.time,final.gravity_values1(:,1),'b'),hold on,
subplot(2,1,2),plot(final.time,final.gravity_values1(:,2),'b'),hold on,

%closer look on April
figure(42)
subplot(2,1,1)
yyaxis left
title('G1-F60')
ylim([-2.9 0.1])
plot(final.time(129601:129600+30*60*24),final.gravity_values1(129601:129600+30*60*24,1),'b')
ylabel('gravity residual [\muGal]')
hold on
subplot(2,1,2)
title('G2-F60')
yyaxis left
plot(final.time(129601:129600+30*60*24),final.gravity_values1(129601:129600+30*60*24,2),'b')
ylabel('gravity residual [\muGal]')
ylim([-5.25 -1.75])
hold on

%plot for resulting timeseries
figure(43)
subplot(2,1,1)
title('G1-F60')
yyaxis left
ylim([-2.9 0.1])
plot(final.time(129601:129600+30*60*24),final.gravity_values1(129601:129600+30*60*24,1),'b')
ylabel('gravity residual [\muGal]')
hold on
subplot(2,1,2)
title('G2-F60')
yyaxis left
plot(final.time(129601:129600+30*60*24),final.gravity_values1(129601:129600+30*60*24,2),'b')
ylabel('gravity residual [\muGal]')
ylim([-5.2 -1.7])
hold on
end

%frequencies to be eliminated
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
%     1/(2*86400)
    %1/(13.1111*86400)
    ]; % day

G1_temp = final.gravity_values1(:,1);
G2_temp = final.gravity_values1(:,2);
%frequency elimination
[x1_hat,G1,Ax1] = FitSin(final.time,G1_temp,f_remaing);
[x2_hat,G2,Ax2] = FitSin(final.time,G2_temp,f_remaing);

if plotting
figure(6)
subplot(2,1,1),plot(final.time,G1 )
title("G1-F60")
ylabel('gravity residual [\muGal]')
legend(["with remaining frequencies","without remaining frequencies"],'Location','southeast')
subplot(2,1,2),plot(final.time,G2 )
title("G2-F60")
ylabel('gravity residual [\muGal]')
legend(["with remaining frequencies","without remaining frequencies"],'Location','southeast')

figure(42)
subplot(2,1,1)
yyaxis right
ylabel('[\muGal]')
ylim([-2.5 0.5])
plot(final.time(129601:129600+30*60*24),Ax1(129601:129600+30*60*24))
legend(["gravity series before frequency elimination","fitted oscillations"],'Location','southeast')
hold off
subplot(2,1,2)
yyaxis right
ylabel('[\muGal]')
ylim([-2 1])
plot(final.time(129601:129600+30*60*24),Ax2(129601:129600+30*60*24))
legend(["gravity series before frequency elimination","fitted oscillations"],'Location','southeast')
hold off

figure(43)
subplot(2,1,1)
yyaxis right
ylabel('[\muGal]')
ylim([-2.7 0.3])
plot(final.time(129601:129600+30*60*24),G1(129601:129600+30*60*24))
legend(["before frequency elimination","after frequency elimination"],'Location','northwest')
hold off
subplot(2,1,2)
yyaxis right
ylabel('[\muGal]')
ylim([-5 -1.5])
plot(final.time(129601:129600+30*60*24),G2(129601:129600+30*60*24))
legend(["before frequency elimination","after frequency elimination"],'Location','southeast')
hold off

%compare psds
[P,f]=plomb(G1,t);
figure(44)
subplot(2,1,2)
semilogy(f,P)
axis([0 5 1e-7 inf])
title('PSD after eliminating frequencies')
xlabel('frequency [cpd]')
ylabel('power spectral density')
end

%final.time_temp = final.time(~idx);
final.time_temp=final.time;
final.gravity_values_temp = [G1 G2];

%% Remove drift
%plot drift
figure(45)
subplot(2,1,1)
final.gravity_values_temp(:,1) = final.gravity_values_temp(:,1) - mean(final.gravity_values_temp(:,1));
final.gravity_values_temp(:,2) = final.gravity_values_temp(:,2) - mean(final.gravity_values_temp(:,2));
plot(final.time,final.gravity_values_temp)
ylabel('gravity residuals [\muGal]')
legend('G1-F60','G2-F60')
title('gravity series before removing drift')

%remove drift
diffG1G2 = final.gravity_values_temp(:,1) - final.gravity_values_temp(:,2);
p = polyfit(seconds(final.time_temp-final.time_temp(1)),diffG1G2,1);
drift = seconds(final.time_temp-final.time_temp(1))*p(1);
final.gravity_values_temp(:,2) = final.gravity_values_temp(:,2) + drift;

figure(45)
subplot(2,1,2)
final.gravity_values_temp(:,1) = final.gravity_values_temp(:,1) - mean(final.gravity_values_temp(:,1));
final.gravity_values_temp(:,2) = final.gravity_values_temp(:,2) - mean(final.gravity_values_temp(:,2));
plot(final.time,final.gravity_values_temp)
title('gravity series after removing drift')
ylabel('gravity residuals [\muGal]')
legend('G1-F60','G2-F60')

%differences between sensors
figure(46)
plot(final.time, final.gravity_values_temp(:,1)-final.gravity_values_temp(:,2))
ylabel('gravity differences [\muGal]')
ylim([-0.5 0.5])

figure(47)
plot(final.time(129601:129600+30*60*24), final.gravity_values_temp(129601:129600+30*60*24,1)-final.gravity_values_temp(129601:129600+30*60*24,2))
yyaxis left
ylabel('gravity differences [\muGal]')
ylim([-0.03 0.09])
yyaxis right
plot(final.time(129601:129600+30*60*24), data.Gravity_cal1.G1_F60(129601:129600+30*60*24), 'Linewidth', 2)
ylabel('gravity values [\muGal]')
ylim([-60 240])

%differences with corrected value
figure(48)
plot(final.time(129601:129600+30*60*24), final.gravity_values_temp(129601:129600+30*60*24,1)-...
    final.gravity_values_temp(129601:129600+30*60*24,2)+...
    0.1/220*40.11*data.Raw.SG.G1_F60(129601:129600+30*60*24))
ylabel('gravity differences [\muGal]')

%hydrology comparison on April 8th
figure(49)
plot(final.time(139681:141120), final.gravity_values_temp(139681:141120,:))
ylabel('gravity differences [\muGal]')
legend('G1-F60','G2-F60')

%hydrology comparison in April
SMSolidus1 = loadSM("BFO_data\Soil_Moisture\SMT100DataSolidus1.csv");
figure(50)
subplot(2,1,1)
% plot(data.Raw.Soil_Moisture.SMT100DataSolidus1.time ,SMSolidus1{1, 3}  )
SM = interp1(seconds(data.Raw.Soil_Moisture.SMT100DataSolidus1.time -...
    data.Raw.Soil_Moisture.SMT100DataSolidus1.time(1)),SMSolidus1{1, 3},seconds(final.time_temp-final.time_temp(1)));
plot(final.time(129601:129600+30*60*24), SM(129601:129600+30*60*24,:))
ylabel('Soil Moisture [%]')
%title("SMT100DataSolidus1")
subplot(2,1,2)
plot(final.time(129601:129600+30*60*24), final.gravity_values_temp(129601:129600+30*60*24,:))
ylabel('gravity residuals [\muGal]')
legend('G1-F60','G2-F60')


if later
%% Final
idxxx = setdiff(final.time,final.time_temp);
final.time_temp = [final.time_temp;idxxx];
final.gravity_values_temp =  [final.gravity_values_temp;[NaN*ones(length(idxxx),1) NaN*ones(length(idxxx),1)]];
[~,idxxxx] = sort(final.time_temp);

final.time = final.time_temp(idxxxx);
final.gravity_values = final.gravity_values_temp(idxxxx,:);

figure(7)
plot(final.time,final.gravity_values(:,1) - final.gravity_values(1,1),'color',[0.8500 0.3250 0.0980],'LineWidth',2	)
hold on
plot(final.time,final.gravity_values(:,2) - final.gravity_values(1,2),'color',[0.4660 0.6740 0.1880],'LineWidth',2)
legend('G1-F60','G2-F60')
%xlabel('time [d]')
ylabel('gravity values [\muGal]')
title('Gavity values in \muGal after correcting the influence of atmospheric pressure and tides (without remaining frequencies)')
ylim([-4 8])
grid on

%% GRACE
%computing means for comparison with GRACE
[monthly_means, daily_means] = GravityResMeans(final.gravity_values,data.Raw.SG.Br2_F60);

run GRACE_data.m

%compare GRACE values with SG
figure(12)
scatter(Grseries(205:214)*1e3-mean(Grseries(205:214)*1e3),monthly_means(1:10),'filled')
xlabel('GRACE gravity residual [\muGal]')
ylabel('average SG gravity residual [\muGal]')
%axis([-2 2.5 -4 4])

%without Love number correction
figure(13)
scatter(Grnoloveseries(205:214)*1e3-mean(Grnoloveseries(205:214)*1e3),monthly_means(1:10),'filled')
xlabel('GRACE gravity residual [\muGal] (without Love number)')
ylabel('average SG gravity residual [\muGal]')
%axis([-2 2.5 -4 4])

%% Soil Moisture
SMSolidus1 = loadSM("BFO_data\Soil_Moisture\SMT100DataSolidus1.csv");
figure(10),
plot(data.Raw.Soil_Moisture.SMT100DataSolidus1.time ,SMSolidus1{1, 3}  )
ylabel('Soil Moisture [%]')
title("SMT100DataSolidus1")

SM = interp1(seconds(data.Raw.Soil_Moisture.SMT100DataSolidus1.time -...
    data.Raw.Soil_Moisture.SMT100DataSolidus1.time(1)),SMSolidus1{1, 3},seconds(final.time_temp-final.time_temp(1)));

%Compare Gravity residuals and SM
figure(14)
%scatter(SM(SM>0),final.gravity_values(SM>0),2)
SMnum=SM(~isnan(final.gravity_values(:,1)));
Gravnum=final.gravity_values(~isnan(final.gravity_values(:,1)),1);
dscatter(SMnum(SMnum>0),Gravnum(SMnum>0));
xlabel('soil moisture [%]')
ylabel('gravity residual [\muGal')
axis tight

figure(15)
yyaxis left
plot(final.time,final.gravity_values(:,1))
ylabel('gravity residual [\muGal]')
yyaxis right
plot(final.time,SM)
ylabel('soil moisture [%]')

%% Pre
len = length(data.Raw.Precipitation.pre)/60;
for i = 1:len-1
    pre_hour(i,1) = sum(data.Raw.Precipitation.pre((i-1)*60+1:(i)*60));
    time_hour(i,1) = data.Raw.Precipitation.time((i-1)*60+1);
end
figure
plot(time_hour,pre_hour,'LineWidth',2,'color','b')
ylabel('Precipitation in mm/hour')
end
