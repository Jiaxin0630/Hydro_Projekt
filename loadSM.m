% put string of data name in

function [E] = loadSM(s)


addpath BFO_data

fileID = fopen(s,'r');
formatSpec = '%s %f %f %f %f';
C = textscan(fileID,formatSpec,'HeaderLines',1,'Delimiter',',');
E=C;


%probe 0
started = 0;
lost = 0;
for i=2:length(C{1,2})
    if started == 2
        if (C{1,2}(i) == 0 || C{1,2}(i)-C{1,2}(i-1)<(-0.3*min(lost+1,20)) || C{1,2}(i)-C{1,2}(i-1)>0.3*min(lost+1,20))
            lost = lost +1;
        elseif lost > 0
            for j = 1:lost
                E{1,2}(i-j)=C{1,2}(i)+(j/(lost+1))*(C{1,2}(i-lost-1)-C{1,2}(i));
            end
            lost = 0;
        end
    elseif started == 0
        if C{1,2}(i) > 0
            started = 1;
        end
    elseif started == 1
        if C{1,2}(i) > 0
            started = 2;
        else
            started = 0;
            E{1,2}(i-1)= 0;
        end
    end
    if (i == length(C{1,2}) && lost > 0)
        for j = 1:lost
            E{1,2}(i-j+1) = C{1,2}(i-lost);
        end
    end
end
% figure()
% plot(E{1,2})
% xlabel('epoche')
% axis tight
% ylabel('soil moisture [%]')
daysofmonth = 31 * ones(12,1);
daysofmonth(2) = 28;
daysofmonth(4:2:6) = 30;
daysofmonth(9:2:11) = 30;

%time
for i=1:length(C{1,2})
    s=C{1,1}{i};
    month=str2num(s(4:5));
    date = month - 1 + (str2num(s(1:2))-1)/daysofmonth(month) + ...
        str2num(s(12:13))/24/daysofmonth(month) + ...
        str2num(s(15:16))/60/24/daysofmonth(month);
    E{1,6}(i)=date;
end
figure()
plot(E{1,6},E{1,2})
%xlabel('epoche')
axis tight
xticklabels({'March','May','July','September','November'})
ylabel('soil moisture [%]')

%probe 1
started = 0;
lost = 0;
for i=2:length(C{1,3})
    if started == 2
        if (C{1,3}(i) == 0 || C{1,3}(i)-C{1,3}(i-1)<(-0.3*min(lost+1,20)) || C{1,3}(i)-C{1,3}(i-1)>0.3*min(lost+1,20))
            lost = lost +1;
        elseif lost > 0
            for j = 1:lost
                E{1,3}(i-j)=C{1,3}(i)+(j/(lost+1))*(C{1,3}(i-lost-1)-C{1,3}(i));
            end
            lost = 0;
        end
    elseif started == 0
        if C{1,3}(i) > 0
            started = 1;
        end
    elseif started == 1
        if C{1,3}(i) > 0
            started = 2;
        else
            started = 0;
            E{1,3}(i-1)= 0;
        end
    end
    if (i == length(C{1,3}) && lost > 0)
        for j = 1:lost
            E{1,3}(i-j+1) = C{1,3}(i-lost);
        end
    end
end
figure()
plot(E{1,6},E{1,3})
%xlabel('epoche')
axis tight
xticklabels({'March','May','July','September','November'})
ylabel('soil moisture [%]')

% %Plotting for April
% figure()
% subplot(2,1,1)
% time = zeros(25372-17771+1,1);
% for i = 17771:25372
%     s=C{1,1}{i};
%     time(i-17770) = str2num(s(1:2)) + (str2num(s(12:13)) + ...
%         str2num(s(15:16))/60)/24;
% end
% plot(time,E{1,2}(17771:25372))
% hold on
% plot(time,E{1,3}(17771:25372))
% hold off
% ylabel('soil moisture [%]')
% %xlabel('time [d]')
% axis([1 31 0 inf])
% %xticks([1 8 15 22 29])
% legend({'probe 0', 'probe 1'})
end