clear;clc;

Rad_dir='';
Albedo_dir='';
ClassType_dir='';
GPP_dir='';
LAI_dir='';
outpath='';

fileFolder=fullfile(GPP_dir);
dir_input=dir(fullfile(fileFolder,'*.tif'));
file_Names={dir_input.name}';
image_num=length(file_Names);

S0=118.1088;   %0.08202MJ/m2.min;1367W/m2;118.1088 MJ/m2.day
omiga_m=[0.6;0.8;0.6;0.8;0.7;0.8;0.8;0.8;0.8;0.9;0.7;0.9;0;0.9;0;0;0];   % clumping index
%1.ENF、2.EBF、3.DNF、4.DBF、5.MF、6.CSH、7.OSH、8.WSA、9.SAV、10.GRA、11.WET、12.CRO、13.Urb、14.CROandNaturalVegetation、15.Snow、16.Barren、17.Water
Emsh_m=[2.04;1.7;1.2;3.07;2.81;7.94;3.36;0.61;0.36;2.98;0;1.08;0;0;0;0;0];
Emsu_m=[1.75;1.34;2.27;1.2;1.1;0.44;3.01;3.26;3.07;2.27;0;3.26;0;0;0;0;0];

for i=1:image_num
    fprintf(num2str(i));
    FileName=file_Names{i};
    
    yearM=str2num(FileName(1:end-4));
    imagedate=Month2Day(yearM);
    Year_num=floor(imagedate/1000);
    Jday_num=mod(imagedate,1000);
    
    LAI_info = geotiffinfo(strcat([LAI_dir,FileName]));
    [LAI,laigeo] = geotiffread(strcat([LAI_dir,FileName])) ;
    LAI(LAI>10 | LAI<0)=0;
    height = LAI_info.Height;
    width = LAI_info.Width;
    
    [ClassType_image,~] = geotiffread(strcat([ClassType_dir,'20100101.tif'])) ;
    [GPP_image,~] = geotiffread(strcat([GPP_dir,FileName])) ;
    GPP_image(GPP_image>3000 | GPP_image<0)=0;
    
    [Albedo_image,Ageo] = geotiffread(strcat([Albedo_dir,num2str(imagedate),'.tif']));
    Albedo_image(Albedo_image>33 | Albedo_image<0)=0;
    
    
    [Rad_image,~] = geotiffread([Rad_dir,'MERRA2_300.tavg1_2d_rad_Nx.',num2str(imagedate),'.tif']);
    Rad_image=imresize(Rad_image,[3600 7200],'bilinear');
    Rad_image(isnan(Rad_image)| Rad_image<0)=0;
    
    Rad_Sg=Rad_image.*0.0864;
    PAR=Rad_image.*0.0864.*0.5;
    
    
    ClassType_image(isnan(ClassType_image)==1)=255;
    
    omiga=zeros(size(ClassType_image));
    Emsh=zeros(size(ClassType_image));
    Emsu=zeros(size(ClassType_image));
    cosza=zeros(size(ClassType_image));
    
    for k=1:size(ClassType_image,1)
        for j=1:size(ClassType_image,2)
            
            ClassID=round(ClassType_image(k,j));
            if ClassID>13
                continue
            else
                omiga(k,j)=omiga_m(ClassID);
                Emsh(k,j)=Emsh_m(ClassID);
                Emsu(k,j)=Emsu_m(ClassID);
            end
        end
    end
    latmax=90;
    lonmix=-180;
    [lon,lat]=lonlat(height,width,latmax,lonmix);
    
    cosza=cosza_cal(lon,lat,Year_num,Jday_num,10,30);
    
    LAIsu=2.*cosza.*(1-exp(-0.5.*omiga.*(LAI./cosza)));
    LAIsh=LAI-LAIsu;
    
    R=Rad_Sg./(S0*8.*cosza);
    
    fR=0.943+0.734.*R-4.9.*R.^2+1.796.*R.^3+2.058.*R.^4;
    fR(R>0.8)=0.13;
    
    fLAI=exp(-0.5.*omiga.*LAI./(0.537+0.025.*LAI));
    
    fCpar=0.07.*omiga.*(1.1-0.1.*LAI).*exp(-cosza);
    fcos=cos(pi/3)./cosza;
    
    Rparsh=((fR.*(1-fLAI))./LAI+((1-fR).*fCpar)).*LAIsh.*(1-Albedo_image);
    Rparsu=((1-fR).*fcos+(fR.*(1-fLAI))./LAI+(1-fR).*fCpar).*LAIsu.*(1-Albedo_image);
    
    GPP_correctpre=GPP_image.*(Emsh.*Rparsh+Emsu.*Rparsu);
    
    
    
    ind=find(isnan(GPP_correctpre) | GPP_correctpre==0 );
    [row,col]=ind2sub(size(GPP_correctpre),ind);
    for index=1:size(row,1)
        GPP_correctpre(sub2ind(size(GPP_correctpre),row(index), col(index)))=GPP_image(row(index), col(index));
    end
    
    fprintf(' %d %s\n',i);
    GeoRef = georasterref('RasterSize', [3600,7200], 'LatitudeLimits', double([-90, 90]), 'LongitudeLimits', double([-180, 180]), 'ColumnsStartFrom', 'north');
    geotiffwrite([outpath,num2str(imagedate),'.tif'],GPP_correctpre,GeoRef);
end