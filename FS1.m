% Script to make Figure S1.  Adapted from Page et al., [2016].
% Page, M. T., Van Der Elst, N., Hardebeck, J., Felzer, K., & Michael, A. J. (2016). Three ingredients for improved global aftershock forecasts: Tectonic region, time‐dependent catalog incompleteness, and intersequence variability. Bulletin of the Seismological Society of America, 106(5), 2290-2301, doi: 10.1785/0120160073.
clear;

% LOAD GLOBAL NEIC CATALOG - code theifed from Morgan Page & slightly modified.
% and strec regionalization (based on Flinn-Engdahl zones and slab data)
load('+MP_Git/NEIC_Catalog_1990-2015.mat')

% FIND MAINSHOCKS
minMainshockMag=6; maxMainshockMag=Inf; 
exclusionDistance=3; excludeDistanceFormat=1;  % 3 fault lengths
exclusionTimeBefore=90; exclusionTimeAfter=10; % in days
maxDepth=50; excludeEarlyCatalog=1; excludeLateCatalog=1;
mainshockIndices = MP_Git.FindMainshockIndices(catalog,minMainshockMag,maxMainshockMag,maxDepth,exclusionDistance,exclusionTimeBefore,exclusionTimeAfter,excludeDistanceFormat,excludeEarlyCatalog,excludeLateCatalog);
regionNames=regionNames(mainshockIndices);

% FIND AFTERSHOCKS OF MAINSHOCKS
startTime=0; endTime=10; % in days
numFaultLengths=3; minDist=5; maxDepth=50;
minMag=4.5; maxMagDiff=Inf;
% 11th column of assignedCatalog gives mainshock index if it's an aftershock, -1 if it's a mainshock, 0 otherwise
assignedCatalog = MP_Git.SortIntoSequences(catalog,mainshockIndices,startTime,endTime,numFaultLengths,minDist,maxDepth,minMag,maxMagDiff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recreate Figure 5 in Page et al. Global Omori parameter paper.
% Find mean maximum likelihood  Omori parameters to STACKED sequences
% within FLINN-ENGDAHL regions with TIME-VARYING Mc.
% Uses an inequality constraint to keep fit higher that fit using all the data above Mcat

G=0.25; Mcat=4.5;  % Parameters for time-varying Mc
b=1;
c0=0.018293; % This is the c-value from the fit to all regions

d=datenum(catalog(:,1),catalog(:,2),catalog(:,3),catalog(:,4),catalog(:,5),catalog(:,6))';
mag=catalog(:,10)';
u=unique(regionNames);
clear exitflag;

u={'ANSR-ABSLDEC'};

for feregion = 1:length(u) % Loop through tectonic regions

  feregionMainshockInd=mainshockIndices(find(strcmp(u(feregion),regionNames)));
  
  aftershockTimes = []; allAftershockTimes = [];
  mainshockMags = [];
  for mainshockIndex=feregionMainshockInd'
    aftInd = find(assignedCatalog(:,11)==mainshockIndex);
    t = d(aftInd)-d(mainshockIndex);
    Mc = max(mag(mainshockIndex)/2 - G -log10(t), Mcat); % time-varying Mc at times of recorded aftershocks
    aboveMcInd = find(mag(aftInd)>=Mc);
    aftershockTimes = [aftershockTimes t(aboveMcInd)]; % times of aftershocks above time-varying Mc
    allAftershockTimes = [allAftershockTimes t]; % times of aftershocks above Mcat
    mainshockMags = [mainshockMags mag(mainshockIndex)];
  end
  
  % ML fit to all earthquakes above Mcat
  Mequiv=(1/b)*log10(sum(10.^(b*mainshockMags)));  % Equivalent mainshock magnitude for stacked sequence
  % Find parameters using all data down to catalog Mc
  LogL = @(x) MP_Git.ComputeLogLikelihood(allAftershockTimes,x(1),x(2),b,c0,Mequiv,Mcat,startTime,endTime);
  x=fminsearch(@(x) -LogL(x), [-2.3 1.2]);
  a0=x(1); p0=x(2); % a and p for inequality constraint
  
  % ML fit with time-varying Mc
  % Log likelihood with inequality constraint
  LogL = @(x) MP_Git.ComputeLogLikelihoodInequality(aftershockTimes,x(1),x(2),b,c0,mainshockMags,Mcat,G,startTime,endTime,a0,p0);
  x=fminsearch(@(x) -LogL(x), [a0 p0]);
  a=x(1); p=x(2);  
  
  % Compute my fits.
  t=0.001:0.001:endTime;
  %[params]=EQ_Rate_Decay_Fit(aftershockTimes,type_flag,[],[]);
  [params_o]=EQ_Rate_Decay_Fit(allAftershockTimes,'Omori',[],[]); [PDFo,CDFo,No]=EQ_Rate_Decay(t,'Omori',params_o);
  [params_e]=EQ_Rate_Decay_Fit(allAftershockTimes,'Exponential',[],[]); [PDFe,CDFe,Ne]=EQ_Rate_Decay(t,'Exponential',params_e);
  [params_s]=EQ_Rate_Decay_Fit(allAftershockTimes,'Stretched',[],[]); [PDFs,CDFs,Ns]=EQ_Rate_Decay(t,'Stretched',params_s);
  [params_c]=EQ_Rate_Decay_Fit(allAftershockTimes,'Cut-off',[],[]); [PDFc,CDFc,Nc]=EQ_Rate_Decay(t,'Cut-off',params_c);
  [params_g]=EQ_Rate_Decay_Fit(allAftershockTimes,'Gamma',[],[]); [PDFg,CDFg,Ng]=EQ_Rate_Decay(t,'Gamma',params_g);
  
  % PLOT 'em up
  aftershockTimes=sort(aftershockTimes); allAftershockTimes = sort(allAftershockTimes);
  figure(feregion); clf
  loglog(allAftershockTimes(2:end),1./diff(allAftershockTimes),'o','Color',[0.95 0.62 0.12], 'DisplayName','AS sample')
  hold on
  loglog(aftershockTimes(2:end),1./diff(aftershockTimes),'bo','LineWidth',2, 'DisplayName','AS sample')
  
  plot(t,10^(a+b*(Mequiv-minMag))*(t+c0).^-p,'k-','LineWidth',2, 'DisplayName','MP Omori fit 1')
  plot(t,10^(a0+b*(Mequiv-minMag))*(t+c0).^-p0,'k--','LineWidth',2, 'DisplayName','MP Omori fit 2')
  plot(t,PDFo,'-', 'Color','#0000FF','LineWidth',1, 'DisplayName','Omori Fit')
  plot(t,PDFe,'-', 'Color','#FF0000','LineWidth',1, 'DisplayName','Exponential Fit')
  plot(t,PDFs,'-', 'Color','#EDB120','LineWidth',1, 'DisplayName','Stretched Fit')
  plot(t,PDFc,'-', 'Color','#FF00FF','LineWidth',1, 'DisplayName','Cut-off Fit')
  plot(t,PDFg,'-', 'Color','#77AC30','LineWidth',1, 'DisplayName','Gamma Fit')
  
  set(gca,'FontSize',16); axis([10^-3 endTime 10^-1 10^7]);
  xlabel('Time since Mainshock (days) ');
  ylabel('Daily aftershock rate ');
  title(strcat(u(feregion), ': a = ', num2str(a), ', p = ',  num2str(p), ', c = ', num2str(c0)))
  legend('Location','southwest')

end





