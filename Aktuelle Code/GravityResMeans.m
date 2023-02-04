% put double with both gravity residual timeseries in

function [monthly_means, daily_means] = GravityResMeans(grav,press)
monthly_means=zeros(12,2);
daily_means=zeros(354,2);
%saving the number of days of each month for averaging over the months
daysofmonth = 31 * ones(12,1);
daysofmonth(2) = 28;
daysofmonth(4:2:6) = 30;
daysofmonth(9:2:11) = 30;
daysofmonth(12)=20;
%set values without pressure readings to NaN
%a=find(data.Raw.SG.Br2_F60==0);
grav(press==0,:)=NaN;
%average over months
n=0;
for i=1:12
    monthly_means(i,:)=nanmean(grav(n+1:n+daysofmonth(i)*60*24,:));
    n=n+daysofmonth(i)*60*24;
end
figure()
plot(monthly_means(:,1)-mean(monthly_means(:,1)))
hold on
plot(monthly_means(:,2)-mean(monthly_means(:,2)))
ylabel('monthly mean gravity residual [\muGal]')
legend({'G1-F60','G2-F60'})
xlabel('time [month]')
axis tight
title('monthly means')
hold off


%daily averaging
n=0;
for i=1:354
%     if(min(isnan(grav(n+1:n+60*24)))>0)
%         daily_means(i,:)=[NaN NaN];
%     else
%         daily_means(i,:)=nanmean(grav(n+1:n+60*24,:));
%     end
    daily_means(i,:)=nanmean(grav(n+1:n+60*24,:));
    n=n+60*24;
end
figure()
plot(daily_means(:,1)-nanmean(daily_means(:,1)))
hold on
plot(daily_means(:,2)-nanmean(daily_means(:,2)))
hold off
ylabel('daily mean gravity residual [\muGal]')
legend({'G1-F60','G2-F60'})
xlabel('time [day]')
axis tight
end