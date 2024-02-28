%%%Projection coordinates and image size must be the same
clear;clc;

Rad_dir='                     \';  %solar incident shortwave radiation data                                 
Albedo_dir='                  \';  %MCD43C3 Albedo product                            
ClassType_dir='               \';  %MCD12Q1 Land Cover Type product                          
GPP_dir='                     \';  %MOD17A2H GPP product                            
LAI_dir='                     \';  %MOD15A2H LAI product
outpath='                     \';  %Correction product
 
fileFolder=fullfile(GPP_dir);
dir_input=dir(fullfile(fileFolder,'*.tif'));
file_Names={dir_input.name}';
image_num=length(file_Names);

%The clumping index,CR,and a for the different vegetation types
%1.ENF、2.EBF、3.DNF、4.DBF、5.MF、6.CSH、7.OSH、8.WSA、9.SAV、10.GRA、11.WET、12.CRO、13.Urb、14.CROandNaturalVegetation、15.Snow、16.Barren、17.Water 
omiga_m=[0.6;0.8;0.6;0.8;0.7;0.8;0.8;0.8;0.8;0.9;0.7;0.9;0;0.9;0;0;0];   % clumping index
a=[0.1;0.19;0.66;0.01;0.01;0.11;0.12;0.01;0.01;0.01;0;0.01;0;0;0;0;0];
Cr=[2.65;1.47;2.08;2.10;2.75;3.61;3.33;3.07;2.70;2.35;0;4.30;0;0;0;0;0];

for i=1:image_num
    FileName=file_Names{i};
    
    yearM=str2num(FileName(1:end-4));
    imagedate=Month2Day(yearM);
    Year_num=floor(imagedate/1000);%Year
    Jday_num=mod(imagedate,1000);  %day
    
    %The LAI 
    LAI_info = geotiffinfo(strcat([LAI_dir,num2str(imagedate),'.tif']));                    
    [LAI,laigeo] = geotiffread(strcat([LAI_dir,num2str(imagedate),'.tif'])) ;
    height = LAI_info.Height;
    width = LAI_info.Width;

    %The MODIS GPP data
    [GPP_MOD,~] = geotiffread(strcat([GPP_dir,FileName])) ;
    
    %The Albedo data
    [Albedo,Ageo] = geotiffread(strcat([Albedo_dir,num2str(imagedate),'.tif']));
    
    %The radiation data
    [Rad_image,~] = geotiffread([Rad_dir,num2str(imagedate),'.tif']);     
    Rad_Sg=Rad_image.*0.0864;      % from W/M2 to MJ/M2/8d  
    PAR=Rad_Sg.*0.5;               % The total incoming photosynthetically active radiation
    
    % The Land Cover Type data
    [ClassType_image,~] = geotiffread(strcat([ClassType_dir,'.tif']));
    omiga=zeros(size(ClassType_image));
    CR=zeros(size(ClassType_image));
    fa=zeros(size(ClassType_image));
    
    for k=1:size(ClassType_image,1)
        for j=1:size(ClassType_image,2)
            ClassID=round(ClassType_image(k,j));
            if ClassID>17                    % Remove non-vegetational cover
                continue
            else
                omiga(k,j)=omiga_m(ClassID); % clumping index  
                CR(k,j)=Cr(ClassID);
                fa(k,j)=a(ClassID);
            end
        end
    end
    latmax=90;
    lonmix=-180;
    [lon,lat]=lonlat(height,width,latmax,lonmix);

    cosza=cosza_cal(lon,lat,Year_num,Jday_num,10,30);%The cosza is the solar zenith angle
    
    % The LAI_su and LAI_sh are the LAI of sunlit and shaded leaves
    LAIsu=2.*cosza.*(1-exp(-0.5.*omiga.*(LAI./cosza)));
    LAIsh=LAI-LAIsu;
    S0=118.1088;            % The S0 is the solar constant (1367 W m-2), 1367 W m-2 = 118.1088 MJ/m2.day 
    R=Rad_Sg./(S0.*cosza);  % The R is the sky clearness index
    %%
    PARdif=PAR.*(0.943+0.734*R-4.9*R.^2+1.796*R.^3+2.058*R.^4);
    PARdif(R>0.8)=PAR(R>0.8)*0.13;
    PARdir=PAR-PARdif;

    PARdifu=PARdif.*exp(-0.5*omiga.*LAI./(0.537+0.025*LAI));

    Cpar=0.07*omiga.*PARdir.*(1.1-0.1*LAI).*exp(-cosza);

    PARsh=(PARdif-PARdifu)./LAI+Cpar;
    PARsu=PARdir.*cos(pi/3)./cosza+PARsh;
    APARsh=(1-Albedo).*PARsh.*LAIsh;
    APARsu=(1-Albedo).*PARsu.*LAIsu;
    
    %%%%% PPFD %%%%
    PPFDsh=((4.55.*PARsh./0.5)./(48*8));
    PPFDsu=((4.55.*PARsu./0.5)./(48*8));%The hour turns to half an hour
    b=1;
    % RPPFDsu and RPPFDsh are the ratio of APAR and the radiation of sunlit and shaded leaves, respectively
    RPPFDsh = ((b./(fa.*PPFDsh+1)).*APARsh)./PAR;
    RPPFDsu = ((b./(fa.*PPFDsu+1)).*APARsu)./PAR;
    %%
    % CTL-MOD17 Model
    GPP_CTL=CR.* GPP_MOD.*(RPPFDsu + RPPFDsh);
    %%
    GeoRef = georasterref('RasterSize', [3600,7200], 'LatitudeLimits', double([-90, 90]), 'LongitudeLimits', double([-180, 180]), 'ColumnsStartFrom', 'north');
    geotiffwrite([outpath,num2str(imagedate),'.tif'],GPP_CTL,GeoRef);
       
    fprintf(' %d %s\n',i);
end

%Calculate the latitude and longitude of each pixel
function [lon,lat]=lonlat(height,width,latmax,lonmix)
    lat=zeros(height,width);
    lon=zeros(height,width);
    %The latitude and longitude of each pixel
    for i = 1 : height
        for j = 1 : width     
            lat(i,j)=latmax-0.05*i;
            lon(i,j)=lonmix+0.05*j;
        end
    end
end
%date conversion
function day =Month2Day(MonthDay)
   day=floor(MonthDay/10000)*1000+floor(datenum(int2str(MonthDay),'yyyymmdd')-datenum(int2str(MonthDay/10000),'yyyy')+1);
end
%Calculate the solar zenith Angle
function [cosza] = cosza_cal(location_time)
%
lat = location_time{1};
lon = location_time{2};

year = location_time{3};
Jday = location_time{4};
hour = location_time{5};
min = location_time{6};

%sun angle
N0 = 79.6764 + 0.2422 * (year-1985) - floor((year-1985) / 4);
t = Jday - N0;
sun_angle = 2 * pi * t / 365.2422;

dec = 0.3723 + 23.2567 * sin(sun_angle) + 0.1149 * sin(2 * sun_angle)-0.1712 * sin(3 * sun_angle)...
    - 0.758 * cos(sun_angle) + 0.3656 * cos(2 * sun_angle) + 0.0201 * cos(3 * sun_angle);
dec = deg2rad(dec);

Eq = 0.0028 - 1.9857 * sin(sun_angle) + 9.9059 * sin(2 * sun_angle) - 7.0924 * cos(sun_angle)...
    - 0.6882 * cos(2 * sun_angle);

timezone = (lon ./ abs(lon)) .* round(abs(lon) ./ 15);

dlon = lon - timezone * 15;
dlon(lon < 0) = timezone * 15 - lon;

local_time = hour + min / 60 + dlon / 15 + Eq / 60;
hour_angle = (local_time - 12) * 15;
hour_angle = deg2rad(hour_angle);

lat = deg2rad(lat);

cosza = sin(lat) .* sin(dec) + cos(lat) .* cos(dec) .* cos(hour_angle);
cosza(cosza < 0) = 0;

end