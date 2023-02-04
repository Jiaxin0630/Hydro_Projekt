function time = timeconvert(time)
time = num2str(time);
time = [str2double(string(time(:,1:4))) str2double(string(time(:,5:6))) str2double(string(time(:,7:8)))];
time = datenum(time(:,1),time(:,2),time(:,3));
end