% Script to plot the real data.  For analysis and QC'ing.
clear;

% Load in a predefined structure for convenience.
tablePData;

% Predefine some values.
GREY=[0.85,0.85,0.85];
ID='SSFS-1993';

% Get the correct index.
i=find(strcmpi({T.ID},ID));
    
% Get the test sample.
load(T(i).file);
cata_flag='main';
if(strcmpi(cata_flag,'sub'))
    cat=S.cat2;
else
    cat=S.cat;
end
T(i).ID


% GR-MFD stats.
[b,~,a,~,~,Mgr,Ngr,ngr]=Bval(cat.mag, T(i).Mc, T(i).dM);
po=[-b,a];
Mgr_fit=[T(i).Mc, max(cat.mag)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Plot the GR-FMD.
figure(1); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
plot(T(i).Mc*[1 1],ylim,'--k');
xlabel('Magnitude (M_L)'); ylabel('Count');
xlim([min(Mgr)-T(i).dM/2 max(Mgr)+T(i).dM/2]); ylim([0.7 1.3*max(Ngr)]);

% Temporal rate plots.
[~,Tbins]=histcounts(cat.time,round(1.5*sqrt(length(cat.time))));
figure(2); clf;
% M-T plot.
ax1=subplot(411);
plot(cat.time,cat.mag,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.time(cat.mag>=T(i).Mc),cat.mag(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Time (Days)'); ylabel('Magnitude (M_L)');
% EQ count rate plot.
ax2=subplot(412);
histogram(cat.time, Tbins, 'FaceColor', GREY); hold on;
histogram(cat.time(cat.mag>=T(i).Mc), Tbins, 'FaceColor', 'g');
xlabel('Time (Days)'); ylabel('Count');
% EQ cumulative count plot.
ax3=subplot(413);
plot(sort(cat.time), 1:length(cat.time),'-k'); hold on;
plot(sort(cat.time(cat.mag>=T(i).Mc)),1:length(cat.time(cat.mag>=T(i).Mc)),'-g' );
xlabel('Time (Days)'); ylabel('Count');
% Volume Rates.
ax4=subplot(414);
if(strcmpi(S.quality,'A'))
    plot(S.inj.time,S.inj.rate);
elseif(strcmpi(S.quality,'B'))
    for j=1:length(S.inj.Ts)
        plot([S.inj.Ts(j) S.inj.Te(j)],[1 1],'-ob'); hold on;
    end
end
xlabel('Time (Days)'); ylabel('Injeciton Rate (m^3/min)');
% Link 'em up.
linkaxes([ax1,ax2],'x');
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax4],'x');

% Spatial plots.
figure(3); clf;
% Map.
ax1=subplot(221);
plot(cat.lon,cat.lat,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.lon(cat.mag>=T(i).Mc),cat.lat(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');
% Latitude cross-section.
ax2=subplot(222);
plot(cat.dep,cat.lat,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.dep(cat.mag>=T(i).Mc),cat.lat(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Depth (km)'); ylabel('Latitude'); zlabel('Depth (km)');
% Longitude cross-section.
ax3=subplot(223);
plot(cat.lon,-cat.dep,'ok','MarkerFaceColor',GREY); hold on;
plot(cat.lon(cat.mag>=T(i).Mc),-cat.dep(cat.mag>=T(i).Mc),'ok','MarkerFaceColor','g');
xlabel('Longitude'); ylabel('Depth (km)');
% 3D plot.
ax4=subplot(224);
plot3(cat.lon,cat.lat,-cat.dep,'ok','MarkerFaceColor',GREY); hold on;
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



    


