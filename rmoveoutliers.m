function [t,timeseries,idx] = rmoveoutliers(t,timeseries)
[~,idx] = rmoutliers(timeseries,'movmedian',hours(5),'SamplePoints',t);
timeseries = timeseries(~idx,:);
t = t(~idx,:);
end