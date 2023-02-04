addpath GRACE_data
addpath SHbundle-master
%addpath uberall-master
%% read data: GRACE-FO
d = dir('GRACE_data\CSR_CSR-Release-06-GFO_96x96_unfiltered\*.gfc');
path = 'GRACE_data\CSR_CSR-Release-06-GFO_96x96_unfiltered\';
for i = 1:length(d)
    filename = [path d(i).name];
    fileID = fopen(filename);
    data_GRACE = textscan(fileID,'%s%f%f%f%f%f%f','headerlines', 35);
    fileID = fopen(filename);
    time = textscan(fileID,...
        '%s%s%s%s',[1,1],'headerlines', 20);
    time_start_temp = time{2};
    time_end_temp = time{4};
    GRACE(i).time_start = time_start_temp;
    GRACE(i).time_end = time_end_temp;
    GRACE(i).data(:,1) = data_GRACE{2};GRACE(i).data(:,2) = data_GRACE{3};
    GRACE(i).data(:,3) = data_GRACE{4};GRACE(i).data(:,4) = data_GRACE{5};
    fclose all;
end

%% Match time
GRACE_ts = char(string({GRACE.time_start}'));
GRACE_te = char(string({GRACE.time_end}'));

GRACE_ts = [str2double(string(GRACE_ts(:,1:4))) str2double(string(GRACE_ts(:,5:6))) ...
    str2double(string(GRACE_ts(:,7:8)))];
GRACE_te = [str2double(string(GRACE_te(:,1:4))) str2double(string(GRACE_te(:,5:6))) ...
    str2double(string(GRACE_te(:,7:8)))];

[GRACE_tyear,GRACE_tmonth] = TimeMatch(GRACE_ts,GRACE_te);


%% Adding degree 1 coefficients
filename = 'GRACE_data\GRACE Technical Note 13.txt';
fileID = fopen(filename);
data_GRACE = textscan(fileID,'%s%f%f%f%f%f%f%f%f','headerlines', 116);
GFZ13b = [data_GRACE{8} data_GRACE{9} data_GRACE{2} data_GRACE{3} data_GRACE{4} data_GRACE{5}];
fclose all;

% plot
idx10 = GFZ13b(:,3)==1 & GFZ13b(:,4)==0;
idx11 = GFZ13b(:,3)==1 & GFZ13b(:,4)==1;

if ~isequal(length(idx10),length(idx11))
    error('the length of the coefficients are not the same !!')
end


C10 = GFZ13b(idx10,5);
C11 = GFZ13b(idx11,5);
S11 = GFZ13b(idx11,6);

t10_start = timeconvert(GFZ13b(idx10,1));
t10_end = timeconvert(GFZ13b(idx10,2));

GFZ13b_ts = datevec(t10_start);
GFZ13b_te = datevec(t10_end);
[GFZ13b_tyear,GFZ13b_tmonth] = TimeMatch(GFZ13b_ts,GFZ13b_te);

GFZ13b_t_C10 = datenum([GFZ13b_tyear GFZ13b_tmonth ones(length(GFZ13b_tyear),1)]);
GFZ13b_t_C11 = GFZ13b_t_C10;


for i = 1:12
    t10_vec = datevec(GFZ13b_t_C10);t11_vec = datevec(GFZ13b_t_C11);
    idx10 = t10_vec(:,2) == i;
    idx11 = t11_vec(:,2) == i;
    mean_C10(i,1) = mean(C10(idx10));
    mean_C11(i,1) = mean(C11(idx11));
    mean_S11(i,1) = mean(S11(idx11));
end

% adding
id = 0;
for i = 2
    idx = find(GFZ13b_tyear == GRACE_tyear(i) & (GFZ13b_tmonth == GRACE_tmonth(i)));

    if isempty(idx)
        id = id+1;
    else
        C10_replace = C10(idx);
        C11_replace = C11(idx);
        S11_replace = S11(idx);
        data_temp = GRACE(i).data;
        GRACE(i).data(data_temp(:,1)==1 & data_temp(:,2)==0,3) = C10_replace;
        GRACE(i).data(data_temp(:,1)==1 & data_temp(:,2)==1,3) = C11_replace;
        GRACE(i).data(data_temp(:,1)==1 & data_temp(:,2)==1,4) = S11_replace;
    end
end


%% C2,0 replacement
filename = 'GRACE_data\NASA GSFC SLR C20 and C30 solutions.txt';
fileID = fopen(filename);
data_GRACE = textscan(fileID,'%f%f%f%f%f%s%s%s%f%f','headerlines', 37);
SLR_C20 = [datenum(datetime(data_GRACE{1},'convertfrom','modifiedjuliandate')) ...
    datenum(datetime(data_GRACE{9},'convertfrom','modifiedjuliandate')) data_GRACE{3}];
fclose all;
SLR_C20_ts = datevec(SLR_C20(:,1));
SLR_C20_te = datevec(SLR_C20(:,2));
[SLR_C20_tyear,SLR_C20_tmonth] = TimeMatch(SLR_C20_ts,SLR_C20_te);

% replacement
id = 0;
for i = 1:length(GRACE_tyear)
    idx = find(SLR_C20_tyear == GRACE_tyear(i) & (SLR_C20_tmonth == GRACE_tmonth(i)));
    if ~isempty(idx)
        C20_replace = SLR_C20(idx,3);
        data_temp = GRACE(i).data;
        GRACE(i).data(data_temp(:,1)==2 & data_temp(:,2)==0,3) = C20_replace;
    else
        id = id+1;

    end
end

for i = 1:length(GRACE)
    GRACE_new(i).year = GRACE_tyear(i);
    GRACE_new(i).month = GRACE_tmonth(i);
    GRACE_new(i).data = GRACE(i).data;
end                                
clearvars GRACE
GRACE = GRACE_new;

%% GRACE data
GM = 6.6743*10^-11*5.9722*10^24;
r = 6378137;
SG_Lat = deg2rad(48.3294);
SG_Lon = deg2rad(8.3284);
g=gaussian(96,1000);
[l,c] = size(GRACE(1).data);
LongMean = zeros(l,c);

for i = 1:length(GRACE)
    LongMean = LongMean + GRACE(i).data;
end
LongMean = mean(LongMean);


for i = 1:length(GRACE)
    clm = GRACE(i).data;
    dclm = clm;
    dclm(:,1:2) = GRACE(i).data(:,1:2);
    dclm_gf = dclm;
    for j = 1:length(dclm_gf)
        dclm_gf(j,3)=dclm_gf(j,3)*g(dclm_gf(j,1)+1);
        dclm_gf(j,4)=dclm_gf(j,4)*g(dclm_gf(j,1)+1);
    end

    sc = clm2sc(dclm_gf);

    Grseries(i)=-gshs_ptw(sc,SG_Lon,SG_Lat,r,r,'quant',...
        'tr','sub_WGS84',false,'GM',GM);
    Grnoloveseries(i)=-gshs_ptwnolove(sc,SG_Lon,SG_Lat,r,r,'quant',...
        'tr','sub_WGS84',false,'GM',GM);
end

figure
t = datetime([GRACE_tyear GRACE_tmonth ones(length(GRACE_tmonth),1)]);
plot(t,(Grseries-mean(Grseries))*1e3,'-o','LineWidth',1)
title('Monthly timeseries of gravity disturbances (minus the average) captured by GRACE Follow-On','FontSize',20)
ylabel('gravity residuals [\muGal]')
ylim([-2 2.5])
%axis([-inf inf -2 2.5])

figure
t = datetime([GRACE_tyear GRACE_tmonth ones(length(GRACE_tmonth),1)]);
plot(t,(Grnoloveseries-mean(Grnoloveseries))*1e3,'-o','LineWidth',1)
title('Monthly timeseries of gravity disturbances (minus the average) captured by GRACE Follow-On','FontSize',20)
ylabel('gravity residuals [\muGal] (without Love numbers)')
ylim([-2 2.5])
%axis([-inf inf -2 2.5])

