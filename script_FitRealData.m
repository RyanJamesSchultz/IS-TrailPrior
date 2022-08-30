% Script to fit real-case data for one case.  Used for analysis and QC'ing.
clear;

% Load in a predefined structure for convenience.
tableTData;

% Predefine some values.
GREY=[0.85,0.85,0.85];
Nf=1e6;
ID='PNR2-g';
lw=1;

% Get the correct index.
i=find(strcmpi({T.ID},ID));
%i=55;
T(i).ID

% Boundary polygons and stuff.
latB=T(i).latB;
lonB=T(i).lonB;
depB=[];
tB=T(i).Tf;
mB=[T(i).Mc Inf];

% Get the test sample.
load(T(i).file);

% Preprocess the data.
[Ts,cat,bounds]=filtCat(S,latB,lonB,depB,tB,mB);
t=linspace(0,bounds(2),Nf); t(1)=[];

% Fit each model to the test sample.
[Sf]=Fit_Trailing(Ts,bounds, Nf,t, T(i).Go,T(i).Ge,T(i).Gs,T(i).Gc,T(i).Gg);
%for j=1:5
%    Sf.Pg(2:end)
%    [Sf]=Fit_Trailing(Ts,bounds, Nf,t, Sf.Po,Sf.Pe,Sf.Ps,Sf.Pc,Sf.Pg);
%end

% Report some values.
fprintf(1,'\n');
fprintf(1,'Waic values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',Sf.Waic(1),Sf.Waic(2),Sf.Waic(3),Sf.Waic(4),Sf.Waic(5));
fprintf(1,'Wbic values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',Sf.Wbic(1),Sf.Wbic(2),Sf.Wbic(3),Sf.Wbic(4),Sf.Wbic(5));
fprintf(1,'log10 KS values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',Sf.KSp(1),Sf.KSp(2),Sf.KSp(3),Sf.KSp(4),Sf.KSp(5));
fprintf(1,'Binned R2 values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',Sf.R2b(1),Sf.R2b(2),Sf.R2b(3),Sf.R2b(4),Sf.R2b(5));
%fprintf(1,'Unbinned R2 values\n');
%fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',Sf.R2r(1),Sf.R2r(2),Sf.R2r(3),Sf.R2r(4),Sf.R2r(5));
fprintf(1,'\nOrdering\n');
fprintf(1,'Omori, Exponential, Stretch, Cut-off, Gamma\n');
n=length(Ts)

% Temporal rate plots.
[~,Tbins]=histcounts(Ts,round(2*sqrt(length(cat.time))));
figure(1); clf;
% M-T plot.
ax1=subplot(411);
plot(S.cat.time-min(tB),S.cat.mag,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.time-min(tB),cat.mag,'ok','MarkerFaceColor','g');
plot(bounds(1)*[1 1], ylim(),'--k');
plot(bounds(2)*[1 1], ylim(),'--k');
xlabel('Days from T-0'); ylabel('Magnitude (M_L)');
% EQ count rate plot.
ax2=subplot(412);
histogram(S.cat.time-min(tB), Tbins, 'FaceColor', GREY); hold on;
histogram(cat.time-min(tB), Tbins, 'FaceColor', 'g');
plot(bounds(1)*[1 1], ylim(),'--k');
plot(bounds(2)*[1 1], ylim(),'--k');
xlabel('Days from T-0'); ylabel('Count');
% EQ cumulative count plot.
ax3=subplot(413);
plot(sort(S.cat.time-min(tB)), 1:length(S.cat.time),'-k'); hold on;
plot(sort(cat.time-min(tB)),1:length(cat.time),'-g' );
plot(bounds(1)*[1 1], ylim(),'--k');
plot(bounds(2)*[1 1], ylim(),'--k');
xlabel('Days from T-0'); ylabel('Count');
% Volume Rates.
ax4=subplot(414);
if(strcmpi(S.quality,'A'))
    plot(S.inj.time-min(tB),S.inj.rate); hold on;
elseif(strcmpi(S.quality,'B'))
    for j=1:length(S.inj.Ts)
        plot([S.inj.Ts(j) S.inj.Te(j)]-min(tB),[1 1],'-ob'); hold on;
    end
end
plot(bounds(1)*[1 1], ylim(),'--k');
plot(bounds(2)*[1 1], ylim(),'--k');
xlabel('Days from T-0'); ylabel('Injeciton Rate (m^3/min)');
% Link 'em up.
linkaxes([ax1,ax2],'x');
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax4],'x');

% Plot trailing fits.
figure(2); clf;
subplot(311);
loglog(Ts(2:end),1./diff(Ts),'ok','DisplayName','Data Sample'); hold on;
loglog(t,Sf.PDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Fit');
loglog(t,Sf.PDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
loglog(t,Sf.PDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
loglog(t,Sf.PDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
loglog(t,Sf.PDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%loglog(bounds(1)*[1 1], ylim(),'--k');
%loglog(bounds(2)*[1 1], ylim(),'--k');
ylabel('Trailing Rate (events/day)'); xlabel('Time (days)');
legend('Location','southwest');
xlim(bounds);
subplot(312);
histogram(Ts,round(2*sqrt(length(Ts))),'Normalization','countdensity','DisplayName','Binned Data Sample'); hold on;
Ylim=ylim();
plot(t,Sf.PDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Fit');
plot(t,Sf.PDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
plot(t,Sf.PDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
plot(t,Sf.PDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
plot(t,Sf.PDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%plot(bounds(1)*[1 1], ylim(),'--k');
%plot(bounds(2)*[1 1], ylim(),'--k');
ylabel('Trailing Rate (events/day)'); xlabel('Time (days)');
legend('Location','northeast');
xlim(bounds); ylim(Ylim*1.1);
subplot(313);
plot(Ts,1:length(Ts),'-k', 'LineWidth',2,'DisplayName','Data Sample'); hold on;
plot(t,Sf.CDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Fit');
plot(t,Sf.CDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
plot(t,Sf.CDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
plot(t,Sf.CDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
plot(t,Sf.CDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%plot(bounds(1)*[1 1], ylim(),'--k');
%plot(bounds(2)*[1 1], ylim(),'--k');
ylabel('Trailing Counts'); xlabel('Time (days)');
legend('Location','southeast');
ylim([0 1.1]*length(Ts));
xlim(bounds);

% Spatial plots.
figure(3); clf;
% Map.
ax1=subplot(221);
plot(S.cat.lon,S.cat.lat,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.lon(cat.mag>=T(i).Mc),cat.lat(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
if(~isempty(latB))
    plot([lonB,lonB(1)],[latB,latB(1)],'-g');
end
xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
% Latitude cross-section.
ax2=subplot(222);
plot(S.cat.dep,S.cat.lat,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.dep(cat.mag>=T(i).Mc),cat.lat(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Depth (km)'); ylabel('Latitude'); zlabel('Depth (km)');
% Longitude cross-section.
ax3=subplot(223);
plot(S.cat.lon,-S.cat.dep,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.lon(cat.mag>=T(i).Mc),-cat.dep(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Longitude'); ylabel('Depth (km)');
% 3D plot.
ax4=subplot(224);
plot3(S.cat.lon,S.cat.lat,-S.cat.dep,'ok','MarkerFaceColor',GREY); hold on;
plot3(cat.lon(cat.mag>=T(i).Mc),cat.lat(cat.mag>=T(i).Mc),-cat.dep(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
%for j=1:length(S.traj)
%    plot3(S(i).traj(j).lon(end),S(i).traj.lat(end),0,'ok','MarkerFaceColor','b');
%    plot3(S(i).traj(j).lon,S(i).traj.lat,-S(i).traj.dep,'-b');
%end
xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
% Link 'em up.
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax1,ax4],'xy');







    


