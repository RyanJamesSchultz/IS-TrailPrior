% Script to make Figure 5.
clear;

% Load in a predefined structure for convenience.
tableTData;

% Predefine some values.
GREY=[0.85,0.85,0.85];
Nf=1e6;
ID='Basel';
lw=1;

% Define the split-times.
%Tc=2:0.5:15;
%Tc=(9:0.1:10)+0.01;
%Tc=[9.77:0.01:9.79,9.81:0.01:9.83];
%Tc=9.81
%Tc=3.5:0.1:5.0;
%Tc=4.0:0.01:4.2;
Tc=4.08;

% Get the correct index.
i=find(strcmpi({T.ID},ID));
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

% Loop over all of the split-times.
for j=1:length(Tc)
    % Fit the peices.
    Tc(j)
    [Sc]=Fit_Trailing2(Ts,bounds, Tc(j), Nf,t, T(i).Go,T(i).Ge,T(i).Gs,T(i).Gc,T(i).Gg);
end

% Plot trailing fits: piecewise and regular.
figure(5); clf;
subplot(121);
loglog(Ts(2:end),1./diff(Ts),'ok','DisplayName','Data Sample'); hold on;
loglog(t,Sf.PDFo,':', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Regular Fit');
loglog(t,Sc.PDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Piecewise Fit');
%loglog(t,Sf.PDFe,':', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Regular Fit');
%loglog(t,Sc.PDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Piecewise Fit');
%loglog(t,Sf.PDFs,':', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Regular Fit');
%loglog(t,Sc.PDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Piecewise Fit');
%loglog(t,Sf.PDFc,':', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Regular Fit');
%loglog(t,Sc.PDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Piecewise Fit');
%loglog(t,Sf.PDFg,':', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Regular Fit');
%loglog(t,Sc.PDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Piecewise Fit');
%loglog(bounds(1)*[1 1], ylim(),'--k');
%loglog(bounds(2)*[1 1], ylim(),'--k');
loglog(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
ylabel('Trailing Rate (events/day)'); xlabel('Time (days)');
legend('Location','southwest');
xlim(bounds);
subplot(122);
semilogx(t,zeros(size(t)),'-k'); hold on;
semilogx(t,residuals1(Ts,t,Sf.CDFo),':','Color','#0000FF','DisplayName','Omori Regular Fit');
semilogx(t,residuals1(Ts,t,Sc.CDFo),'-','Color','#0000FF','DisplayName','Omori Piecewise Fit');
%semilogx(t,residuals1(Ts,t,Sf.CDFe),':','Color','#FF0000','DisplayName','Exponential Regular Fit');
%semilogx(t,residuals1(Ts,t,Sc.CDFe),'-','Color','#FF0000','DisplayName','Exponential Piecewise Fit');
%semilogx(t,residuals1(Ts,t,Sf.CDFs),':','Color','#EDB120','DisplayName','Stretched Regular Fit');
%semilogx(t,residuals1(Ts,t,Sc.CDFs),'-','Color','#EDB120','DisplayName','Stretched Piecewise Fit');
%semilogx(t,residuals1(Ts,t,Sf.CDFc),':','Color','#FF00FF','DisplayName','Cut-off Regular Fit');
%semilogx(t,residuals1(Ts,t,Sc.CDFc),'-','Color','#FF00FF','DisplayName','Cut-off Piecewise Fit');
%semilogx(t,residuals1(Ts,t,Sf.CDFg),':','Color','#00FF00','DisplayName','Gamma Regular Fit');
%semilogx(t,residuals1(Ts,t,Sc.CDFg),'-','Color','#00FF00','DisplayName','Gamma Piecewise Fit');
plot(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
xlabel('Time (days)'); ylabel('Cumulative Density Residual (-)');
xlim([min(t) max(t)]);
legend('Location','southwest');




% Subroutines.
function res = residuals1(xd,xf,yf)
  % Local subroutine that calculates the CDF residuals (on fitted x-axis).
  
  yd=interp1(xd+1e-10*(1:length(xd))',1:length(xd),xf,'nearest','extrap');
  res=(yf-yd)/length(xd);
  
end