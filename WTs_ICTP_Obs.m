% Atmospheric circulation states based on K-means analysis of atm circulations 
% AG Munoz (IRI Columbia U; CMC Universidad del Zulia) - agmunoz@iri.columbia.edu
% Project: IRAP
% First edition: Feb 4, 2013
% Last edition: Nov 24, 2014

%Key questions here: 
%(1) filter the annual cycle? No need to filter it for the present study
%EOF-filtering is applied.
%
% Details:
% The code uses NNRP data from a particular domain (sdomain) to define
% clusters representing circulation states, usually in terms of geopotential heigth
% at 500 mb, but this is modifiable. Composites of moisture fluxes, rainfall,
% and lightning are computed for each cluster. 
%
% Maps correspond to a bigger domain (bdomain) or could be customizable (pdomain). 
%
% Data will be read via DAP and written to a _ICTPObs.mat file. If the domains/variables are
% not changed, the user doesn't need to read the data again from the DAP
% server.
%
% Note (only if using chi): year 2009 has a lot of missing values for chi as for Nov '14. 
%An easy way to solve the problem is just to make yeare=2008.


%%%%%START OF USER-MODIFIABLE SECTION%%%%%%%%%%%%

disp('Start...');
% set working directory
clear all
% set working directory
cd /Users/agms/Documents/Angel/GFDL/WTs/Weather_within_climate
%addpath /usr/local/bin
addpath /Users/agms/Documents/MATLAB/m_map
%parpool('local')

%Read data via OpenDAP?
down=0;   %1=yes; 0=no (this assumes the data is available in _ICTPObs.mat format; *not* NetCDF!)

%Select variable for clustering:
var=1;   %options available are 
         % 1= z500 (geopotential height at 500mb)
         % 2= chi  (velocity potential)

%Define temporal parameters:
seasons='Dec';  %start
midmon ='Jan';  %middle month
seasone='Feb';  %end
yeari=1981; %first year (MUST BE >=1995 AND <2005!)
yeare=2010; %last year  (MUST BE >2006!)    %Note: for chi there're (or used to be) missing values in 2009!!!!

%Define spatial parameters:
%sdomain (for clusters)
% slonmin=263.75;
% slonmax=306.25;
% slatmin=5;
% slatmax=30;
slonmin=-105;
slonmax=-69;
slatmin=30;
slatmax=50;
%bdomain (for plotting)
blonmin=-180;
blonmax=180;
blatmin=0;
blatmax=90;
%pdomain (custom plotting)
pdlatmin=slatmin;
pdlatmax=slatmax;
pdlonmin=slonmin;
pdlonmax=slonmax;

%Define cluster parameters:
minclust=5; %min num of clust, typically 2
maxclust=5; %typically 10
varfract=0.95; %total required variance explained by the EOF pre-filtering in the k-means algorithm
nclust=5; % choice of k (use classifiability index and physics to determine this value)




%%%END OF USER-MODIFIABLE SECTION (DO NOT MODIFY ANYTHING BELOW THIS LINE!!)%%%%%


%% Rainfall 
% Due to the way IRIDL has divided the dataset (retro + present), we need
% to read them in two pieces and add a bridge between them.
% Future version: use appendstream at DL;
%                 sdomain is a subset of bdomain, so manage that in Matlab
%% Rainfal  YEARI-2005 CPC Unified Precipitation RETRO + DJF 2006-YEARE CPC Unified Precipitation Realtime
if down==1
months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};
%RETRO:
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.RETRO/.rain/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%20' num2str(yeari) '-2005%29RANGEEDGES/dods'];
%Xp=double(ncread(iridl,'X'));
%Yp=double(ncread(iridl,'Y'));


pr1 = double(ncread(iridl,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr1));  %get dims
pr1 = squeeze(pr1);
pr1 = permute(pr1,[2 1 3]);  %we want lat first
 

% This is just to plot it if necessary (plot section comes later; this is a
% test)
% figure(1); clf
%  Xmat=repmat(Xp',length(Yp),1); Ymat=repmat(Yp,1,length(Xp));
%  m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);  
%   m_coast('patch',[.7 .7 .7],'edgecolor','none');
%   m_grid;
%   hold on
%   [cs,h]=m_contour(Xmat,Ymat,squeeze(pr1(:,:,1000)));
  
pr1=reshape(pr1,nlatp*nlonp,ndat2)'; 

%Juncture: Dec 2005 + Jan-Feb 2006
iridl1=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.RETRO/.rain/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '%202005%29RANGEEDGES/missing_value/-999./def/dods'];
pr10 = double(ncread(iridl1,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr10));
pr10 = squeeze(pr10);
pr10 = permute(pr10,[2 1 3]);
pr10=reshape(pr10,nlatp*nlonp,ndat2)';

%iridl2='http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.u/T/%28Jan-Feb%202006%29RANGE/X/-65/-53.5/RANGEEDGES/Y/-36/-25/RANGEEDGES/P/%28500%29VALUES/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/T/%28Jan-Feb%202006%29RANGEEDGES/X/-65/-53.5/RANGEEDGES/Y/-36/-25/RANGEEDGES%5BT%5DregridLinear//missing_value/-999./def/dods';
iridl2=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' midmon '-' seasone '%202006%29RANGEEDGES/missing_value/-999./def/dods'];
pr11 = double(ncread(iridl2,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr11));
pr11 = squeeze(pr11);
pr11 = permute(pr11,[2 1 3]);
pr11=reshape(pr11,nlatp*nlonp,ndat2)'; 

%Realtime:
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%202006-' num2str(yeare) '%29RANGEEDGES/missing_value/-999./def/dods'];
pr2 = double(ncread(iridl,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr2));
pr2 = squeeze(pr2);
pr2 = permute(pr2,[2 1 3]);
pr2=reshape(pr2,nlatp*nlonp,ndat2)';

%Concatenate RETRO, juncture and Realtime:
prj= cat(1,pr10,pr11);
pr_SESA = cat(1,pr1,prj,pr2);

clear prj pr1 pr11 pr10 pr12 pr2 iridl1 iridl2
%% Rainfal DJF YEARI-2005 CPC Unified Precipitation RETRO + DJF 2006-YEARE CPC Unified Precipitation Realtime

%RETRO:
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.RETRO/.rain/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%20' num2str(yeari) '-2005%29RANGEEDGES/dods'];
Xp=double(ncread(iridl,'X'));
Yp=double(ncread(iridl,'Y'));


pr1 = double(ncread(iridl,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr1));  %get dims
pr1 = squeeze(pr1);
pr1 = permute(pr1,[2 1 3]);  %we want lat first
  
pr1=reshape(pr1,nlatp*nlonp,ndat2)'; 

%Juncture: Dec 2005 + Jan-Feb 2006
iridl1=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.RETRO/.rain/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '%202005%29RANGEEDGES/missing_value/-999./def/dods'];
pr10 = double(ncread(iridl1,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr10));
pr10 = squeeze(pr10);
pr10 = permute(pr10,[2 1 3]);
pr10=reshape(pr10,nlatp*nlonp,ndat2)';

%iridl2='http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.u/T/%28Jan-Feb%202006%29RANGE/X/-65/-53.5/RANGEEDGES/Y/-36/-25/RANGEEDGES/P/%28500%29VALUES/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/T/%28Jan-Feb%202006%29RANGEEDGES/X/-65/-53.5/RANGEEDGES/Y/-36/-25/RANGEEDGES%5BT%5DregridLinear//missing_value/-999./def/dods';
iridl2=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' midmon '-' seasone '%202006%29RANGEEDGES/missing_value/-999./def/dods'];
pr11 = double(ncread(iridl2,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr11));
pr11 = squeeze(pr11);
pr11 = permute(pr11,[2 1 3]);
pr11=reshape(pr11,nlatp*nlonp,ndat2)'; 

%Realtime:
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP/.CPC/.UNIFIED_PRCP/.GAUGE_BASED/.GLOBAL/.v1p0/.REALTIME/.rain/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/T/2.0/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%202006-' num2str(yeare) '%29RANGEEDGES/missing_value/-999./def/dods'];
pr2 = double(ncread(iridl,'rain'));
[nlonp nlatp ndat2]=size(squeeze(pr2));
pr2 = squeeze(pr2);
pr2 = permute(pr2,[2 1 3]);
pr2=reshape(pr2,nlatp*nlonp,ndat2)';

%Concatenate RETRO, juncture and Realtime:
prj= cat(1,pr10,pr11);
pr = cat(1,pr1,prj,pr2);

clear prj pr1 pr11 pr10 pr12 pr2 iridl1 iridl2
save -v7.3 precip_ICTPObs.mat pr pr_SESA Xp Yp 
else
load precip_ICTPObs.mat pr pr_SESA Xp Yp
end
disp('Rainfall (Unified CPC, 0.5 deg) has been read and storaged ');
%% CIRCULATION VARIABLE
%This is sdmain (for computing clusters):
if down==1

  if var==1
  % z850
  iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.phi/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/P/%28500%29VALUES/dods'];
  phi = double(ncread(iridl,'phi'));
  elseif var==2
  % chi850
  iridl=['http://iridl.ldeo.columbia.edu/expert/%28http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/daily/sigma/chi%29readdods/.chi/%28http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/daily/surface_gauss/pres%29readdods/.pres/lev/mul/lev/85000/toS/time//T/renameGRID/lon//X/renameGRID/lat//Y/renameGRID/T/1.5/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(slonmin) '/' num2str(slonmax) '/RANGEEDGES/Y/' num2str(slatmin) '/' num2str(slatmax) '/RANGEEDGES/T/0/FillNaN/dods'];
  phi = double(ncread(iridl,'chi'));
  %phi ( phi==-32766998 )=0;
  else
  disp('Unkown option for cluster variable');
  end

[nlon2 nlat2 ndat2]=size(squeeze(phi));
phi = squeeze(phi);
phi = permute(phi,[2 1 3]);

%T=double(ncread(iridl,'T'));
X=double(ncread(iridl,'X'));
Y=double(ncread(iridl,'Y'));

phi=reshape(phi,nlat2*nlon2,ndat2)';

save -v7.3 var4clust_ICTPObs.mat phi X Y ndat2 
else
load var4clust_ICTPObs.mat X Y phi ndat2 
end
disp('Circulation variable (NNRPv2, 2.5 deg) has been read and storaged ');
%% SECTION TO COMPUTE CLUSTERS via K-means
rng(1); % initialize random numbers for reproducible results

CI=NaN*ones(maxclust,1); K=NaN*ones(ndat2,maxclust);
for kk=minclust:maxclust
  rng(1); % initialize random numbers for reproducible results
  [CI(kk),K(:,kk)]=kmeans_ci(phi,'s',varfract,kk,100);
  disp(['Calculated k-means for K=' num2str(kk) ' yielding CI of ' num2str(CI(kk))])
end

% %
% % % Plot CI Index
% figure(1); clf
% plot(CI,'k-','Linewidth',1.5)
% hold on 
% set(gca,'FontSize',14)
% xlabel('No. of Clusters'); ylabel('Classifiability Index')
% title(['Classifiability Index - ' num2str(yeari) '-' num2str(yeare) ' ' num2str(varfract) ' var'])
% grid on
% axis([minclust maxclust 0.7 1.02])
save -v7.3 CI_95p_z500_obs_ICTPObs.mat CI K


indx=K(:,nclust);
%% MORE DATA FOR PLOT SECTION
%%%This subsection is in charge of DAPing data to plot bdomain/pdomain
if down==1

  if var==1
  % z500
  iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.phi/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/P/%28500%29VALUES/dods'];
  phi = double(ncread(iridl,'phi'));
  elseif var==2
  % chi850
  iridl=['http://iridl.ldeo.columbia.edu/expert/%28http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/daily/sigma/chi%29readdods/.chi/%28http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/daily/surface_gauss/pres%29readdods/.pres/lev/mul/lev/85000/toS/time//T/renameGRID/lon//X/renameGRID/lat//Y/renameGRID/T/1.5/runningAverage/T//pointwidth/1.0/def/pop/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/T/0/FillNaN/dods'];
  phi = double(ncread(iridl,'chi'));
  %phi ( phi==-32766998 )=0;
  else
  disp('Unkown option for cluster variable');
  end

Xv=double(ncread(iridl,'X'));
Yv=double(ncread(iridl,'Y'));

[nlonv nlatv ndat2]=size(squeeze(phi));
phi = squeeze(phi);
phi = permute(phi,[2 1 3]);
phi=reshape(phi,nlatv*nlonv,ndat2)';
%phi=phi(1:double(lseas)*nseas,:);   

disp('Circulation variable (NNRPv2, 2.5 deg) for bdomain has been read and storaged ');

% moisture fluxes
% uq 
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.qa/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.u/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/mul/P/low/high/RANGE/P/500/1000/definite-integral/dods'];
int_dP = double(ncread(iridl,'int_dP'));
[nlon2 nlat2 ndat2]=size(squeeze(int_dP));
int_dP = squeeze(int_dP);
uq = permute(int_dP,[2 1 3]);
uq=reshape(uq,nlat2*nlon2,ndat2)';
T=double(ncread(iridl,'T'));
X=double(ncread(iridl,'X'));
Y=double(ncread(iridl,'Y'));

% vq
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.qa/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.DAILY/.Intrinsic/.PressureLevel/.v/T/%28' seasons '-' seasone '%20' num2str(yeari) '-' num2str(yeare) '%29RANGE/X/' num2str(blonmin) '/' num2str(blonmax) '/RANGEEDGES/Y/' num2str(blatmin) '/' num2str(blatmax) '/RANGEEDGES/mul/P/low/high/RANGE/P/500/1000/definite-integral/dods'];
clear int_dP
int_dP = double(ncread(iridl,'int_dP'));
[nlon2 nlat2 ndat2]=size(squeeze(int_dP));
int_dP = squeeze(int_dP);
vq = permute(int_dP,[2 1 3]);
vq=reshape(vq,nlat2*nlon2,ndat2)';

disp('Moisture fluxes (NNRPv2, 2.5 deg) have been read and storaged ');

save -v7.3 plotvars_ICTPObs.mat phi uq vq nlat2 nlon2 nlatp nlonp nlonv nlatv Xv Yv X Y
else
load plotvars_ICTPObs.mat phi uq vq nlat2 nlon2 nlatp nlonp nlonv nlatv Xv Yv X Y
end
disp('Data have been read and storaged ');


disp('The fun is about to start...');
%
%I don't understand why precip is 1 day shorter
nt=2616; %2797; %ndat2-1;     %in case a subset is needed
ndat2=nt;
phi=phi(1:nt,:);
pr=pr(1:nt,:);
uq=uq(1:nt,:);
vq=vq(1:nt,:);
indx=K(1:nt,nclust);

%%
%Computing anomalies
clear uqcompa vqcompa z3compa prcompa 
for k=1:nclust  
  clear kk
  kk=find(indx==k);
  nday(k)=length(kk);
  %computing anomaly fields:
  uqcompa(k,:)=squeeze(nanmean(uq(kk,:),1) - nanmean(uq,1));
  vqcompa(k,:)=squeeze(nanmean(vq(kk,:),1) - nanmean(vq,1));
  z3compa(k,:)=squeeze(nanmean(phi(kk,:),1) - nanmean(phi,1));
  prcompa(k,:)=squeeze(nanmean(pr(kk,:),1) - nanmean(pr,1));
  prcompo(k,:)=squeeze(nanmean(pr(kk,:),1));
end
save -v7.3 h500_NENA_ICTPObs.mat z3compa nday Xv Yv Xp Yp prcompa prcompo

uqcompa=reshape(uqcompa,nclust,nlat2,nlon2);
vqcompa=reshape(vqcompa,nclust,nlat2,nlon2);
z3compa=reshape(z3compa,nclust,nlatv,nlonv);
prcompa=reshape(prcompa,nclust,nlatp,nlonp);
prcompo=reshape(prcompo,nclust,nlatp,nlonp);
[dumb,kkplot]=sort(nday,'descend');
prclim = reshape(nanmean(pr,1),nlatp,nlonp);
z3clim = reshape(nanmean(phi,1),nlatv,nlonv);    %For Xiaosong's tests

kkplot = [4 5 2 1 3];
%% PLOT SECTION
%%plot composites
disp('Generating plots ');

%Circulation variable and moisture flux anomalies
clear varcompa X2 Y2 Xv2 Yv2
[X2,Y2]=meshgrid(X,Y);
[Xv2,Yv2]=meshgrid(Xv,Yv);
if var==2
    %Let's interpolate chi to the standard NNRP grid
    for k=1:nclust
    varcompa(k,:,:) = interp2(Xv2, Yv2, squeeze(z3compa(k,:,:)),X2,Y2);
    end
else
    varcompa=z3compa;
end
%clear z3compa 
Xmat=repmat(X',length(Y),1); Ymat=repmat(Y,1,length(X));
figure(2); clf
labs={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
set(gca,'FontSize',14)
%ns=1.; %n-sigmas
pct_th=66;  %percentile threshold for plotting arrows
scal=1.5;
colormap((b2r(-10,10)))
pdlatmin=blatmin;
pdlatmax=blatmax;
pdlonmin=blonmin;
pdlonmax=blonmax;
for kplot=1:nclust
  map=squeeze(varcompa(kkplot(kplot),:,:));
  map(isnan(map))=0;
  umap=squeeze(uqcompa(kkplot(kplot),:,:))*scal; %
  vmap=squeeze(vqcompa(kkplot(kplot),:,:))*scal; %

   um  =prctile(reshape(umap,size(umap,1)*size(umap,2),1),pct_th);    
   vm  =prctile(reshape(vmap,size(umap,1)*size(umap,2),1),pct_th);
   
   umap ( umap>-um & umap<um  )=NaN;
   vmap ( vmap>-vm & vmap<vm  )=NaN;
   
   h=subplot(1,nclust,kplot);
   position=get(h,'position');

  %m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]); 
  m_proj('stereographic','longitude',[-90],'latitude',[90],'radius',[80]);
  m_coast('patch',[.7 .7 .7],'edgecolor','none');
%m_grid;
m_grid('xticklabels',[]);

  hold on
  %Add a square to show sdomain
  %m_line(-59,-30.5,'marker','square','markersize',65,'color','black','linewidth',1.5);  hold on

  if var==1 % z850
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',0.9); %caxis([-15 15]);
  hold on
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',1.05); %caxis([-15 15]);
  elseif var==2 %chi
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',0.9);
  else
  disp('Unkown option for cluster variable');
  end
  
  %clabel(cs,h,'fontsize',12);
  hold on
  %m_quiver(Xmat,Ymat,umap,vmap,0,'color','black');


  %m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]); 
  m_proj('stereographic','longitude',[-90],'latitude',[90],'radius',[80]);
  
  %set(gca,'FontSize',12)
  title([labs{kplot} ' WT ' num2str(kplot) '  (' num2str(round(nday(kkplot(kplot))/ndat2*100)) '% of days)'])
end
colormap((b2r(-200,200)))
h=colorbar;
set(h, 'Position', [.92 .235 .02 .69])
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,['WTs - z500 - DJF ' num2str(yeari) '-' num2str(yeare) ' ' num2str(varfract) ' var'],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

   %In what follows, 2*sqrt(um*um+vm*vm) was selected to provide ~100 g/kg
   %m/s. It may differ in other studies
   %For m_vec we have SCALE, LAT, LON, MAG, etc)
  %[hpv5, htv5] = m_vec(1, 1, -40, 2*sqrt(um*um+vm*vm)*scal, 0, 'black', 'key', '100 g kg^{-1} m s^{-1}');
  %set(htv5,'FontSize',12);
  
orient landscape
disp('WTs figures saved')

% %Precipitation anomaly
% pdlatmin=slatmin;
% pdlatmax=slatmax;
% pdlonmin=slonmin;
% pdlonmax=slonmax;
% Xmat=repmat(Xp',length(Yp),1); Ymat=repmat(Yp,1,length(Xp));
% figure(3); clf
% labs={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
% colormap(flipud(b2r(-10,10)))
% for kplot=1:nclust
%   map=squeeze(prcompa(kkplot(kplot),:,:));
%   map(isnan(map))=0;
%   h=subplot(1,nclust,kplot);
%   position=get(h,'position');
%   m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);  
%   m_coast('patch',[.7 .7 .7],'edgecolor','none');
% %m_grid;
% m_grid('xticklabels',[]);
%   %Add a square to show sdomain
%   %m_line(-59,-30.5,'marker','square','markersize',65,'color','black','linewidth',1.5);
%   hold on
%   [cs,h]=m_contour(Xmat,Ymat,map); caxis([-1 1]);
%   m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);
%   title([labs{kplot} ' WT ' num2str(kplot) '  (' num2str(round(nday(kkplot(kplot))/ndat2*100)) '% of days)'])
% end
% h=colorbar;
% set(h, 'Position', [.92 .235 .02 .69])
% %ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% %text(0.5, 1,['Rainfall - WTsz500 - DJF ' num2str(yeari) '-' num2str(yeare) ' ' num2str(varfract) ' var'],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
% 

%%------------------------------------------
%Precipitation anomaly
pdlatmin=slatmin;
pdlatmax=slatmax;
pdlonmin=slonmin;
pdlonmax=slonmax;
Xmat=repmat(Xp',length(Yp),1); Ymat=repmat(Yp,1,length(Xp));
figure(3); clf
labs={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
%colormap(flipud(b2r(-10,10)))
for kplot=1:nclust
  map=squeeze(prcompa(kkplot(kplot),:,:));
  map(isnan(map))=0;
  h=subplot(1,nclust,kplot);
  position=get(h,'position');
  m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);  
 
  %m_grid('xtick',[pdlonmin:5:pdlonmax],'ytick',[pdlatmin:10:pdlatmax]); 
  m_grid('xticklabels',[]);

  m_coast('linewidth',2,'color','black');
        hold on

            m_plus = map;
            m_plus(map<=0) = NaN;
            m_neg = map;
            m_neg(map>0) = NaN;

            [cs,h]=m_contour(Xmat,Ymat,m_neg,[-2:0.5:0],'linewidth',1.5,'color','blue','LineStyle','--'); %caxis([-100 100]);  %orig steo 0.2
            clabel(cs,h,'fontsize',10,'color','black');
            hold on
            [cs,h]=m_contour(Xmat,Ymat,m_plus,[0:0.5:2],'linewidth',1.5,'color','red'); %caxis([-100 100]); %orig steo 0.2
            clabel(cs,h,'fontsize',10,'color','black');

        %following two lines are used to plot stat sig values (**masking**)
%         map=squeeze(fsigp(kplot,:,:)).*map;  
%         map(map==0)=NaN;      %masking
%         %map(isnan(map))=0; 
%         [cs,h]=m_contourf(Xmatp,Ymatp,map,[-100:80:20]);
%         caxis([-100 20]);
%         colormap(gray)
%         hold on

         %Choose one of the following projections:
        %Stereographic (NH)
        %m_proj('stereographic','longitude',[-90],'latitude',[90],'radius',[80]); 
        %m_grid('xticklabels',[],'yticklabels',[]);

  m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);
   title([labs{kplot} ' WT ' num2str(kplot) '  (' num2str(round(nday(kkplot(kplot))/ndat2*100)) '% of days)'])
end
%h=colorbar;


orient landscape
disp('Rainfall figures saved')

Xmat=repmat(X',length(Y),1); Ymat=repmat(Y,1,length(X));
figure(4); clf
labs={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
set(gca,'FontSize',14)
%ns=1.; %n-sigmas
pct_th=66;  %percentile threshold for plotting arrows
scal=0.3;
colormap((b2r(-10,10)))
pdlatmin=15;
pdlatmax=70;
pdlonmin=-120;
pdlonmax=-60;
for kplot=1:nclust
  map=squeeze(varcompa(kkplot(kplot),:,:));
  map(isnan(map))=0;
  umap=squeeze(uqcompa(kkplot(kplot),:,:))*scal; %
  vmap=squeeze(vqcompa(kkplot(kplot),:,:))*scal; %

   um  =prctile(reshape(umap,size(umap,1)*size(umap,2),1),pct_th);    
   vm  =prctile(reshape(vmap,size(umap,1)*size(umap,2),1),pct_th);
   
   umap ( umap>-um & umap<um  )=NaN;
   vmap ( vmap>-vm & vmap<vm  )=NaN;
   
   h=subplot(1,nclust,kplot);
   position=get(h,'position');

  m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]); 
  m_coast('patch',[.7 .7 .7],'edgecolor','none');
  %m_grid;
  m_grid('xticklabels',[]);
  
  hold on
  %Add a square to show sdomain
  %m_line(-59,-30.5,'marker','square','markersize',65,'color','black','linewidth',1.5);  hold on

  if var==1 % z850
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',0.9); %caxis([-15 15]);
  hold on
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',1.05); %caxis([-15 15]);
  elseif var==2 %chi
  [cs,h]=m_contour(Xmat,Ymat,map,'linewidth',0.9);
  else
  disp('Unkown option for cluster variable');
  end
  
  clabel(cs,h,'fontsize',12);
  hold on
  m_quiver(Xmat,Ymat,umap,vmap,0,'color','black');


  m_proj('Equidistant Cylindrical','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]); 
  
  %set(gca,'FontSize',12)
  title([labs{kplot} ' WT ' num2str(kplot) '  (' num2str(round(nday(kkplot(kplot))/ndat2*100)) '% of days)'])
end
colormap((b2r(-200,200)))
h=colorbar;
set(h, 'Position', [.92 .235 .02 .69])
%ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

%text(0.5, 1,['WTs - z500 - DJF ' num2str(yeari) '-' num2str(yeare) ' ' num2str(varfract) ' var'],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

   %In what follows, 2*sqrt(um*um+vm*vm) was selected to provide ~100 g/kg
   %m/s. It may differ in other studies
   %For m_vec we have SCALE, LAT, LON, MAG, etc)
  %[hpv5, htv5] = m_vec(1, 1, -40, 2*sqrt(um*um+vm*vm)*scal, 0, 'black', 'key', '100 g kg^{-1} m s^{-1}');
  %set(htv5,'FontSize',12);
  
orient landscape
disp('WTs figures saved')

%Precipitation climatology
pdlatmin=slatmin;
pdlatmax=slatmax;
pdlonmin=slonmin;
pdlonmax=slonmax;
Xmat=repmat(Xp',length(Yp),1); Ymat=repmat(Yp,1,length(Xp));
figure(5); clf
colormap(jet)
  map=squeeze(prclim);
  %map(isnan(map))=0;
  m_proj('lambert','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);  
  m_coast('color',[0. 0. 0.],'linewidth',2.5);
  m_grid;
  %Add a square to show sdomain
  %m_line(-59,-30.5,'marker','square','markersize',65,'color','black','linewidth',1.5);
  hold on
  [cs,h]=m_contourf(Xmat,Ymat,map); %caxis([-1 1]);
  m_proj('lambert','longitudes',[pdlonmin pdlonmax],'latitudes',[pdlatmin pdlatmax]);
h=colorbar;
set(gca,'FontSize',14)
  title('A) Observations')

%set(h, 'Position', [.92 .195 .02 .69])
%set(h, 'Position', [.92 .235 .02 .69])
%ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

save -v7.3 precip_obs_ICTPObs.mat prclim Xp Yp 
orient landscape
disp('Rainfall clim figures saved')

save -v7.3 precons_obs_ICTPObs.mat prcompo Xp Yp 

%%
%Others (normally, don't modify this part)
nseas=yeare-yeari;   %number of seasons (1980 is Dic1979+JanFeb1980)
lseas=90; %int8(ndat2/nseas);   %length of the season (DJF in our case); lseas MUST BE EVEN INTEGER, adjust accordingly!!!!  2790  31

%%
%%WTs stats
%indx=K(1:double(lseas)*nseas,nclust); 
indx=K(1:lseas*nseas,nclust);    
clear nday ndayperseas
indx=reshape(indx,lseas,nseas);
for iyr=1:nseas
  for k=1:nclust  
  clear kk kk1 kk2
            kk=find(indx(:,iyr)==k);
            ndayperseas(iyr,k)=length(kk);
  end
end
%Re-ordering:
for k1=1:nclust
      ndayperseaso(:,k1) = ndayperseas(:,kkplot(k1));
end

figure(6); clf
%plot(ndayperseas,'Linewidth',1.5)
bar(ndayperseaso,'stack')
%bar(ndayperseas)
legend('1','2','3','4','5','6','7','8')
set(gca,'FontSize',14)
xlabel('DJF'); ylabel('Frequency (days)')
xlim([0.5 nseas+.5])
ylim([0 lseas+0.5])
set(gca,'XTick',1:5:nseas)
set(gca,'XTickLabel',{'1981','1986','1991','1996','2001','2006','2011','2016'})
%title('Occupation Frequency by Year')
grid off
cmap=colormap;

%Klee diagram!!!
indx=K(1:lseas*nseas,nclust); 
indx=reshape(squeeze(indx(1:nseas*lseas,:)),lseas,nseas);
 clear indxo
% %Re-ordering:
for nd=1:lseas
    for ns=1:nseas
        for k=1:nclust
            if indx(nd,ns) == kkplot(k) 
                indxo(nd,ns)= k;
            end
        end
    end
end

figure(7); clf
colormap((b2r(-5,5)))
%bar3(indx(:,4:5))
imagesc(indxo(:,:))
%axis([1 90 1 2 ])
%set(gca,'XTick',0.5:1:2.5)
%set(gca,'XTickLabel',{'82-83','83-84'})
set(gca,'FontSize',15)
xlabel('DJF'); ylabel('Calendar Day')
set(gca,'XTick',1:5:nseas)
set(gca,'XTickLabel',{'1981','1986','1991','1996','2001','2006','2011','2016'})
%colormap(coll);%colorbar
set(gca,'FontSize',15)
%title('Klee Diagram - Threshold=1 mm ')

%% occ freq by cal day
clear ndayperCalDay F xx Fstd Fo algo upper_limit lower_limit
for iday=1:lseas   
  for k=1:nclust  
  clear kk
  kk=find(indx(iday,:)==k);    %full period
  ndayperCalDay(iday,k)=length(kk);
  end
end

fn=11;  
a = 1;
b = ones(1,fn)/fn;
F = filter(b,a,ndayperCalDay); 
F=cat(1,NaN*ones((fn-1)/2,nclust),F(fn:end,:),NaN*ones((fn-1)/2,nclust));
% Fstd=nanstd(F,0,1);
% 
% % Approximative
% for k =1:nclust
%   upper_limit(:,k)=F(:,k)+0.83*Fstd(k);
%   lower_limit(:,k)=F(:,k)-0.83*Fstd(k);
% end
% %F=filter0(ndayperCalDay,[inf 20]);

figure(77); clf
xx=1:90;
%col=colormap('default');
coll=colormap(cmap);

%Re-ordering:
for k1=1:nclust
      Fo(:,k1) = F(:,kkplot(k1));
end
colormap((b2r(-5,5)))
bar(Fo,'stack')
legend('1','2','3','4','5','6','7','8')
set(gca,'FontSize',14)
xlabel('Calendar Day'); ylabel('Frequency (years)')
xlim([6 84])
ylim([0 nseas])
%set(gca,'XTick',1:5:nseas)
%set(gca,'XTickLabel',{'1981','1986','1991','1996','2001','2006','2011','2016'})
%title('Occupation Frequency by Year')
grid off


%Transition Probabilities
indx=K(1:ndat2,nclust);    
clear nday 
for k=1:nclust  
  clear kk
  kk=find(indx==k);
  nday(k)=length(kk);
end
%State transition probabilities
%chunk into 2 day sequences
temp1=reshape(K(1:ndat2,nclust),2,ndat2/2);    % cols are 2-day seqs
temp2=reshape(K(2:ndat2-1,nclust),2,(ndat2-2)/2); % same, staggered by 1 day

for k1=1:nclust
    for k2=1:nclust
      tr(k1,k2) = length(find(temp2(1,:)==k1 & temp2(2,:)==k2));  ...
    end
end

for k1=1:nclust
    for k2=1:nclust
      ntran(k1,k2) = length(find(temp1(1,:)==k1 & temp1(2,:)==k2)) ...
                + length(find(temp2(1,:)==k1 & temp2(2,:)==k2)); % k1 to k2
    end
end

colsum=sum(ntran,1);
rowsum=sum(ntran,2);

for k1=1:nclust
    for k2=1:nclust
      ptran(k1,k2) = round(100*ntran(k1,k2)/rowsum(k1));
    end
end
%Re-ordering:
for k1=1:nclust
    for k2=1:nclust
      ptrano(k1,k2) = ptran(kkplot(k1),kkplot(k2));
    end
end
coll=colormap;

figure(8); clf
hb=bar3(ptrano)
colormap(jet)
for i = 1:length(hb)
     zdata = get(hb(i),'Zdata');
     set(hb(i),'Cdata',zdata)
end
axis([0.5 5.5 0.5 5.5 0 90])
set(gca,'FontSize',15)
title('Transition Matrix')



save -v7.3 stats_obs_ICTPObs.mat ptrano ndayperseaso indxo Fo
%% SEASONAL SOURCES OF PREDICTABILITY:



%% SST composites for DJF when WT freq exceeds p-th percentile
clear sst sstsig Xsst Ysst sstcomp sstcompa ssty
iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.NOAA/.NCDC/.ERSST/.version4/.sst/T/%28' num2str(yeari-1) '-' num2str(yeare) '%29VALUES/T/%28' seasons '-' seasone '%29seasonalAverage/dods'];
ncid = netcdf.open ( iridl );
 %[name,xtype,dimids,natts] = netcdf.inqVar(ncid,1);
 %name
sst = double(ncread(iridl,'sst'));;
Xsst= double(ncread(iridl,'X'));
Ysst=double(ncread(iridl,'Y'));
%
p=80;  %percentile!!
thresh = prctile(ndayperseas,p);

%for k=6:6%nclust
for k=1:nclust
    iyears{k} = find(ndayperseas(:,k)>thresh(k));
    sstcomp(:,:,k)=mean(sst(:,:,1,iyears{k}),4);
    sstcompa(:,:,k)=mean(sst(:,:,1,iyears{k}),4) - mean(sst(:,:,1,:),4);
end


% %%------------------------------------------
% %UNCOMMENT IF T-TEST IS DESIRED
% sig=0.05;  %confidence limit
% for k=1:nclust
%     iyears{k} = find(ndayperseas(:,k)>thresh(k));
%      %Simple t-test for statistical significance of anomalies (comparing to
%      %null anomaly)
%      ssty=squeeze(squeeze(sst(:,:,1,iyears{k})));  %SST!!!, not anomalies.
%      for ilon=1:length(Xsst)
%         for ilat=1:length(Ysst)
%             sstsig(ilon,ilat,k)=1-ttest(ssty(ilon,ilat,:)-mean(sst(ilon,ilat,1,:),4),0.,'Alpha',sig); 
%         end
%      end
% end
% %%%------------------------------------------


%%%------------------------------------------
%%UNCOMMENT IF HYPERGEOMETRIC TEST IS DESIRED
%Hypergeometric test of statistical significance of anomalies
%Ingredients:
%Zeroth: confidence limit
sig=0.05;
%First: population size
%dime=size(T);
pop=ndat2;%dime(1);
%Second: number of events for each gridbox in POPULATION (event: anom>0,
%anom<0, etc): Kevp_p,n (p for positive, n for negative)
clear Kevp_p Kevp_n
for ilon=1:length(Xsst)
       for ilat=1:length(Ysst)
           Kevp_p(ilon,ilat)=nansum(squeeze(squeeze(sst(ilon,ilat,1,:)))-mean(sst(ilon,ilat,1,:),4)>0);  %counting positive anom in POP
           Kevp_n(ilon,ilat)=nansum(squeeze(squeeze(sst(ilon,ilat,1,:)))-mean(sst(ilon,ilat,1,:),4)<0);  %counting negative anom in POP
       end
end
%Third: sample size ( this is just: length(iyears{k}), k=1..nclust ).
%AND
%Fourth: number of events for each gridbox in SAMPLE (event: anom>0, anom<0,
%etc): : Kevs_p,n (p for positive, n for negative)
clear sample Kevs_p Kevs_n sstsig_p sstsig_n sstsig
for k=1:nclust
    iyears{k} = find(ndayperseas(:,k)>thresh(k));
    sample(k) = length(iyears{k});
    ssty=squeeze(squeeze(sst(:,:,1,iyears{k})));  %SST!!!, not anomalies.
    for ilon=1:length(Xsst)
       for ilat=1:length(Ysst)
           Kevs_p(ilon,ilat,k)=nansum(ssty(ilon,ilat,:)-mean(sst(ilon,ilat,1,:),4)>0);  %counting positive anom in SAMP
           Kevs_n(ilon,ilat,k)=nansum(ssty(ilon,ilat,:)-mean(sst(ilon,ilat,1,:),4)<0);  %counting negative anom in SAMP
           %Hypergeometric!
           % probability of getting at least Kevs_p,n occurrences of the phenomenon --by 
           % chance-- out of a set of n events  P(X>=Kevs|pop,Kevp,sampl):
           sstsig_p(ilon,ilat,k)=1-hygecdf(Kevs_p(ilon,ilat,k)-1,pop,Kevp_p(ilon,ilat),sample(k));   %positives
           sstsig_n(ilon,ilat,k)=1-hygecdf(Kevs_n(ilon,ilat,k)-1,pop,Kevp_n(ilon,ilat),sample(k));   %negatives
           %Proceed to mask
           %positives:
            if(sstsig_p(ilon,ilat,k)<=sig)
                sstsig_p(ilon,ilat,k)=1;
            else
                sstsig_p(ilon,ilat,k)=NaN;
            end
            %negatives:
            if(sstsig_n(ilon,ilat,k)<=sig)
                sstsig_n(ilon,ilat,k)=1;
            else
                sstsig_n(ilon,ilat,k)=NaN;
            end
            %final:
            if((sstsig_p(ilon,ilat,k)+sstsig_n(ilon,ilat,k))==2)
                sstsig(ilon,ilat,k)=1;
            else
                sstsig(ilon,ilat,k)=NaN;
            end
       end
    end
end
clear sstsig_n sstsig_p
%%%End of Hypergeom test
%%------------------------------------------


%So sstsig=1 means statistically significant anomalies.
%%%PLOTS

Xmat=repmat(Xsst',length(Ysst),1)'; Ymat=repmat(Ysst,1,length(Xsst))';
figure(10); clf
colormap((b2r(-0.5,0.5)))
for kplot=1:nclust
  map=squeeze(sstcompa(:,:,kkplot(kplot)));
  %following two lines are used to plot stat sig values
  map=squeeze((sstsig(:,:,kkplot(kplot))).*map);  
  map(map==0)=NaN;      %masking
  subplot(1,nclust,kplot)
  m_proj('Equidistant Cylindrical','longitudes',[120 360],'latitudes',[-30 70]);
  m_coast('patch',[.7 .7 .7],'edgecolor','none');
  m_grid;
  hold on
  [cs,h]=m_contourf(Xmat,Ymat,map,[-1:0.1:1]); caxis([-.5 .5]); %colorbar
% clabel(cs,h,'manual');
  title(['WT ' num2str(kplot) '  (' num2str(length(iyears{kkplot(kplot)})) ' years)' ])
end
h=colorbar;
set(h, 'Position', [.92 .42 .02 .2])

%%
%%Indices
%ENSO - Ni�o3.4
another=['http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.nino/.EXTENDED/.NINO34/T/%28' num2str(yeari-1) '-' num2str(yeare) '%29VALUES/T/%28' seasons '-' seasone '%29seasonalAverage/dods'];
ncid = netcdf.open ( another );
%[name,xtype,dimids,natts] = netcdf.inqVar(ncid,1);
%name
NINO34 = netcdf.getVar(ncid,1);
%n34=reshape(NINO34,3,nseas);
n34m=NINO34';
%n34m=squeeze(mean(n34));
clear ENSOcorr PVAL ENSOcorro
for k=1:nclust
    [ENSOcorr(k), PVAL(k)]=corr(n34m',ndayperseas(1:nseas,k));
    %[ENSOcorr(k), PVAL(k)]=corr(n34m(end-60:end)',ndayperseas(end-60:end,k))
end

%% bootstrapping - scramble Nino34 timeseries
nsamp=10000; clear signif signifo
for i=1:nsamp
 ransample=n34m(randperm(nseas));
 for k=1:nclust
   ENSOcorrBoot(k,i)=corr(ransample',ndayperseas(1:nseas,k));
 end
end
for k=1:nclust
    ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',2.5)); 
    ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
    ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
    ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',97.5)); 
%     ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',10)); 
%     ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
%     ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
%     ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',90)); 

 if(ENSOcorr(k)<ENSOcorrBoot_5(k) || ENSOcorr(k)>ENSOcorrBoot_95(k))
      signif{k}='*';
  else
      signif{k}=' ';
  end
  if(ENSOcorr(k)<ENSOcorrBoot_2p5(k) || ENSOcorr(k)>ENSOcorrBoot_97p5(k))
      signif{k}='**';
  end
end

figure(11); clf
subplot(1,3,1);
%reordering:
for kk=1:nclust
    ENSOcorro(kk)=ENSOcorr(kkplot(kk));
end
bar(ENSOcorro,'FaceColor', [0.5 0.5 0.5])
for kk=1:nclust
    text(kk-0.3,0.45,signif{kkplot(kk)},'FontSize',18)
end

set(gca,'FontSize',14)
xlabel('WT'); ylabel('Ni�o3.4 Correlation' )
axis([0 nclust+1 -.5 .5])
set(gca,'YTick',[-.5:.1:.5])
%title('Correlations between Nino3.4 (NDJ) and WT Freq. ')
%grid on
%orient landscape
%print -dpdf CorrelENSO_OBS.pdf

%
%PNA
another=['http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.CPC_Indices/.NHTI/.PNA/T/%28' num2str(yeari-1) '-' num2str(yeare) '%29VALUES/T/%28' seasons '-' seasone '%29seasonalAverage/dods'];
ncid = netcdf.open ( another );
%[name,xtype,dimids,natts] = netcdf.inqVar(ncid,1);
%name
NINO34 = netcdf.getVar(ncid,1);
%n34=reshape(NINO34,3,nseas);
n34m=NINO34';
%n34m=squeeze(mean(n34));
clear ENSOcorr PVAL ENSOcorro
for k=1:nclust
    [ENSOcorr(k), PVAL(k)]=corr(n34m',ndayperseas(1:nseas,k));
    %[ENSOcorr(k), PVAL(k)]=corr(n34m(end-60:end)',ndayperseas(end-60:end,k))
end

%% bootstrapping - scramble timeseries
nsamp=10000; clear signif signifo
for i=1:nsamp
 ransample=n34m(randperm(nseas));
 for k=1:nclust
   ENSOcorrBoot(k,i)=corr(ransample',ndayperseas(1:nseas,k));
 end
end
for k=1:nclust
    ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',2.5)); 
    ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
    ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
    ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',97.5)); 
%     ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',10)); 
%     ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
%     ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
%     ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',90)); 

 if(ENSOcorr(k)<ENSOcorrBoot_5(k) || ENSOcorr(k)>ENSOcorrBoot_95(k))
      signif{k}='*';
  else
      signif{k}=' ';
  end
  if(ENSOcorr(k)<ENSOcorrBoot_2p5(k) || ENSOcorr(k)>ENSOcorrBoot_97p5(k))
      signif{k}='**';
  end
end


subplot(1,3,2);
%reordering:
for kk=1:nclust
    ENSOcorro(kk)=ENSOcorr(kkplot(kk));
end
bar(ENSOcorro,'FaceColor', [0.5 0.5 0.5])
for kk=1:nclust
    text(kk-0.3,0.45,signif{kkplot(kk)},'FontSize',18)
end

set(gca,'FontSize',14)
xlabel('WT'); ylabel('PNA Correlation' )
axis([0 nclust+1 -.5 .5])
set(gca,'YTick',[-.5:.1:.5])
%title('Correlations between Nino3.4 (NDJ) and WT Freq. ')
%grid on
%orient landscape
%print -dpdf CorrelPNA_OBS.pdf

%
%NAO
another=['http://iridl.ldeo.columbia.edu/SOURCES/.Indices/.CPC_Indices/.NHTI/.NAO/T/%28' num2str(yeari-1) '-' num2str(yeare) '%29VALUES/T/%28' seasons '-' seasone '%29seasonalAverage/dods'];
ncid = netcdf.open ( another );
%[name,xtype,dimids,natts] = netcdf.inqVar(ncid,1);
%name
NINO34 = netcdf.getVar(ncid,1);
%n34=reshape(NINO34,3,nseas);
n34m=NINO34';
%n34m=squeeze(mean(n34));
clear ENSOcorr PVAL ENSOcorro
for k=1:nclust
    [ENSOcorr(k), PVAL(k)]=corr(n34m',ndayperseas(1:nseas,k));
    %[ENSOcorr(k), PVAL(k)]=corr(n34m(end-60:end)',ndayperseas(end-60:end,k))
end

%% bootstrapping - scramble timeseries
nsamp=10000; clear signif signifo
for i=1:nsamp
 ransample=n34m(randperm(nseas));
 for k=1:nclust
   ENSOcorrBoot(k,i)=corr(ransample',ndayperseas(1:nseas,k));
 end
end
for k=1:nclust
    ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',2.5)); 
    ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
    ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
    ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',97.5)); 
%     ENSOcorrBoot_2p5(k)=(prctile(ENSOcorrBoot(k,:)',10)); 
%     ENSOcorrBoot_5(k)=(prctile(ENSOcorrBoot(k,:)',5)); 
%     ENSOcorrBoot_95(k)=(prctile(ENSOcorrBoot(k,:)',95)); 
%     ENSOcorrBoot_97p5(k)=(prctile(ENSOcorrBoot(k,:)',90)); 

 if(ENSOcorr(k)<ENSOcorrBoot_5(k) || ENSOcorr(k)>ENSOcorrBoot_95(k))
      signif{k}='*';
  else
      signif{k}=' ';
  end
  if(ENSOcorr(k)<ENSOcorrBoot_2p5(k) || ENSOcorr(k)>ENSOcorrBoot_97p5(k))
      signif{k}='**';
  end
end


subplot(1,3,3);
%reordering:
for kk=1:nclust
    ENSOcorro(kk)=ENSOcorr(kkplot(kk));
end
bar(ENSOcorro,'FaceColor', [0.5 0.5 0.5])
for kk=1:nclust
    text(kk-0.3,0.45,signif{kkplot(kk)},'FontSize',18)
end
set(gca,'FontSize',14)
xlabel('WT'); ylabel('NAO Correlation' )
axis([0 nclust+1 -.5 .5])
set(gca,'YTick',[-.5:.1:.5])
%grid on
%orient landscape
%print -dpdf CorrelNAO_OBS.pdf


%%
% SUB-SEASONAL 
%% MJO stratification
% lag-composites of RMM index (daily, 1974/6/1 - 2011/10/3)
clear indx
rmm=dlmread('MJO_RMM12.txt', ' ');
rmm=rmm(:,1:7);
% select DJF 1981--2010 (29 yrs)
rmmDJF=rmm((rmm(:,2)>11 | rmm(:,2)<=2),:);
rmmDJF=rmmDJF(633:3249,:);
rmmDJF_phase=rmmDJF(1:nseas*lseas,6);
% select last 15 days of Nov:
%%%% Feb 14 - 28 in non-leap yrs
%%%% Feb 15 - 29 in leap yrs
rmmF=NaN*ones(nseas,15);
for yr=1981:2009
    rmmF(yr-1980,:)=rmm(rmm(:,1)==yr & rmm(:,2)==11 & rmm(:,3)>15,6)';
%   if(mod((yr-1972),4)>0) % non-leap yrs (Munich Olympics as key)
%       %disp([num2str(yr) 'non-leap'])
%     rmmF(yr-1978,:)=rmm(rmm(:,1)==yr & rmm(:,2)==2 & rmm(:,3)>=14,6)';
%   else % leap years
%       %disp([num2str(yr) 'LEAP'])
%     rmmF(yr-1978,:)=rmm(rmm(:,1)==yr & rmm(:,2)==2 & rmm(:,3)>=15,6)';
%   end
end

% Cassou(2008, Nature)
% Table of contingency between the MJO phases (rows) and the 
% North Atlantic weather regimes (columns). For each MJO phase, 
% I plot the anomalous percentage occurrence of a given regime
% as a function of lag in days (with regimes lagging MJO phases).
% The 0% value means that the MJO phase is not discriminative for
% the regime whose occurrence is climatological. A 100% value would
% mean that this regime occurs twice as frequently as its climatological
% mean; -100% means no occurrence of this regime. 

%%
 indx=K(1:nseas*lseas,nclust); %

% marginals
for k=1:nclust  
  ndayWT(k)=length(find(indx==k));
end
for mjop=1:8  
  ndayMJOP(mjop)=length(find(rmmDJF_phase==mjop));
end
% contingency between MJO phase cat & WT
clear cgcy
for k=1:nclust  % (row)
  for mjop=1:8  % (col)
      cgcy(k,mjop)=length(find(rmmDJF_phase==mjop & indx==k));
  end
end

clear cgcyboot cgcya
% bootstrap, by resampling indx
nsamp=500;
for i=1:nsamp
 scrambledi=indx(randperm(length(indx)));
 for k=1:nclust  % (row)
  for mjop=1:8  % (col)
      cgcyboot(k,mjop,i)=length(find(rmmDJF_phase==mjop & scrambledi==k));
  end
 end  
end
cgcyboot_5=round(prctile(cgcyboot,5,3));   %orig 5!
cgcyboot_50=round(prctile(cgcyboot,50,3));
cgcyboot_95=round(prctile(cgcyboot,95,3));  %orig 95!
cgcya=cgcy-squeeze(cgcyboot_50);
cgcya5=cgcyboot_5-squeeze(cgcyboot_50);
cgcya95=cgcyboot_95-squeeze(cgcyboot_50);

% figure(69);clf
% colormap((b2r(-10,10)))
% cgcya(cgcya>=cgcya5 & cgcya<=cgcya95)=0;
% imagesc(cgcya,[-40 40]); 
% %axis([0 nclust+1 -.5 .5])
% set(gca,'FontSize',14)
% xlabel('MJO Phase');
% ylabel('WT')
% title('All years');


% compute at lag -- with MJO leading the WT index -- up to 2 weeks
% WT index - DJF
% MJO index - extends back into Nov
rmmDJF_phase=rmmDJF_phase(1:nseas*lseas,:);
indx=indx(1:nseas*lseas,:);
rmmDJF_phase=reshape(rmmDJF_phase,lseas,nseas);
rmm_ext=cat(1,rmmF',rmmDJF_phase);
cgcyLead=NaN*ones(8,nclust,15);
for lead=1:15 % (incl zero lag)
  % first make a suitably lagged MJO index
  rmm_ext1=rmm_ext(17-lead:17-lead+(lseas-1),:); % for lead 1, start at first day of Dec    
  rmm_ext1=reshape(rmm_ext1,lseas*nseas,1);
  for k=1:nclust  % (row)
    for mjop=1:8  % (col)
      cgcyLead(k,mjop,lead)=length(find(rmm_ext1==mjop & indx==k));
    end
  end
end

% bootstrap lagged WTs, by resampling indx
nsamp=500;
cgcyLead_boot5=zeros(nclust,8,15);
cgcyLead_boot50=zeros(nclust,8,15);
cgcyLead_boot95=zeros(nclust,8,15);

for lead=1:15 % (incl zero lag)
  rmm_ext1=rmm_ext(17-lead:17-lead+(lseas-1),:); % for lead 1, start at first day of Dec    
  rmm_ext1=reshape(rmm_ext1,lseas*nseas,1);
  for k=1:nclust  % (row)
    for mjop=1:8  % (col)
      for i=1:nsamp
        scrambledi=indx(randperm(length(indx)));
        bootsamp(i)=length(find(rmm_ext1==mjop & scrambledi==k));
      end
      cgcyLead_boot5(k,mjop,lead)=prctile(bootsamp,5); 
      cgcyLead_boot50(k,mjop,lead)=prctile(bootsamp,50); 
      cgcyLead_boot95(k,mjop,lead)=prctile(bootsamp,95);
    end
  end
end

%% plot matrix
figure(12); clf
colormap((b2r(-20,20)))
for k=1:nclust
subplot(1,nclust,k)
datamap=squeeze(cgcyLead(kkplot(k),:,:))'-squeeze(cgcyLead_boot50(kkplot(k),:,:))';
c5map=squeeze(cgcyLead_boot5(kkplot(k),:,:))'-squeeze(cgcyLead_boot50(kkplot(k),:,:))';
c95map=squeeze(cgcyLead_boot95(kkplot(k),:,:))'-squeeze(cgcyLead_boot50(kkplot(k),:,:))';
datamap(datamap>=c5map & datamap<=c95map)=NaN;

imagesc(datamap,[-20 20]); cmap=colormap;
%imagesc(datamap); cmap=colormap;
cmap(1,:)=[1 1 1]; colormap(cmap);
%colorbar
grid on
%set(gca,'FontSize',14)
xlabel('MJO Phase')
if(k==1) ylabel('MJO Lead Time (days)'), end
set(gca,'YTick',[1 3 6 8 11 15])
set(gca,'YTickLabel',{'0','2','5','7','10','14'})
title(['WT ' num2str(k)])
end
%subplot(2,3,6)
h=colorbar;
set(h, 'Position', [.92 .16 .02 .74])
orient landscape
set(gcf,'Renderer','zbuffer')