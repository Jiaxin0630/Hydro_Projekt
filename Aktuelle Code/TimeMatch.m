% function [y,m] = TimeMatch(s,e)
% year_s  = s(:,1);
% year_e = e(:,1);
% month_s  = s(:,2);
% month_e = e(:,2);
% day_s = s(:,3);
% day_e = e(:,3);
% month_before = 999;
% for i = 1:length(year_s)
%     if month_s(i) == month_e(i)
%         y(i,1) = year_s(i);
%         m(i,1) = month_e(i);
%     elseif abs(day_s(i)-30) > abs(day_e(i)-1) && month_s(i) ~= month_before
%         y(i,1) = year_s(i);
%         m(i,1) = month_s(i);
%     else
%         y(i,1) = year_e(i);
%         m(i,1) = month_e(i);
%     end
%     month_before = m(i,1);
% end
% end


function [y,m] = TimeMatch(s,e)
year_s  = s(:,1);
year_e = e(:,1);
month_s  = s(:,2);
month_e = e(:,2);
day_s = s(:,3);
day_e = e(:,3);
month_before = 999;
for i = 1:length(year_s)
    if isequal(month_s(i),month_e(i)) && (day_e(i) - day_s(i) >= 10)
        y(i,1) = year_s(i);
        m(i,1) = month_e(i);
    elseif ~isequal(month_s(i),month_e(i)) && ~isequal(month_s(i),month_before)
        y(i,1) = year_s(i);
        m(i,1) = month_s(i);
    elseif ~isequal(month_s(i),month_e(i)) && isequal(month_s(i),month_before) && (month_e(i) - month_s(i) == 1)
        y(i,1) = year_e(i);
        m(i,1) = month_e(i);
    elseif ~isequal(month_s(i),month_e(i)) && isequal(month_s(i),month_before) && (month_e(i) - month_s(i) ~= 1)
        if month_s == 12
            y(i,1) = year_s(i);
        else
            y(i,1) = year_e(i);
        end
        m(i,1) = month_s(i)+1;
    else
        error("check the matching !!")
    end
    month_before = m(i,1);
    
end
if ~any(diff(m))
    fprintf('check the data !!!\n')
end
end