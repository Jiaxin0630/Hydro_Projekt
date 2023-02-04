function [t,timeseries2] = rmoveoutliers(t,timeseries,timeseries2)
[~,idx] = rmoutliers(timeseries,'movmedian',hours(5),'SamplePoints',t);
idxx = find(idx == 1);
i=1;
while i<=length(idxx)
    j=1;
    while true
        %how many more outliers are succeding
        if(i+j <= length(idxx) && idxx(i+j)==idxx(i)+j)
            j=j+1;
        else
            %interpolate linearly for all of them
            for k = 0:j-1
                timeseries2(idxx(i)+k,:)=timeseries2(idxx(i)-1,:)+(k+1)/j*...
                    (timeseries2(idxx(i)+j,:)-timeseries2(idxx(i)-1,:));
            end
            break;
        end
    end
    %jump to next cluster
    i=i+j;
end
%     vq1 = interp1(t([idxx(i)-1,idxx(i)+1]), timeseries2([idxx(i)-1,idxx(i)+1],1),t(idxx(i)));
%     vq2 = interp1(t([idxx(i)-1,idxx(i)+1]), timeseries2([idxx(i)-1,idxx(i)+1],2),t(idxx(i)));
%     timeseries2(idxx(i),1) = vq1;
%     timeseries2(idxx(i),2) = vq2;
end         
