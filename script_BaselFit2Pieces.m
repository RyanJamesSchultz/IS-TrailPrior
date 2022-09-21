% Script to fit real-case data for one case.  Used for analysis and QC'ing.
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

% Define vectors for all of the performance metrics.
ll_o=zeros(size(Tc)); ll_e=ll_o; ll_s=ll_o; ll_c=ll_o; ll_g=ll_o;
aic_o=ll_o; aic_e=ll_o; aic_s=ll_o; aic_c=ll_o; aic_g=ll_o;
bic_o=ll_o; bic_e=ll_o; bic_s=ll_o; bic_c=ll_o; bic_g=ll_o;
ksp_o=ll_o; ksp_e=ll_o; ksp_s=ll_o; ksp_c=ll_o; ksp_g=ll_o;

% Loop over all of the split-times.
for j=1:length(Tc)
    % Fit the peices.
    Tc(j)
    [Sc]=Fit_Trailing2(Ts,bounds, Tc(j), Nf,t, T(i).Go,T(i).Ge,T(i).Gs,T(i).Gc,T(i).Gg);
    ll_o(j)=Sc.LL(1);    ll_e(j)=Sc.LL(2);   ll_s(j)=Sc.LL(3);   ll_c(j)=Sc.LL(4);   ll_g(j)=Sc.LL(5);
    aic_o(j)=Sc.AIC(1); aic_e(j)=Sc.AIC(2); aic_s(j)=Sc.AIC(3); aic_c(j)=Sc.AIC(4); aic_g(j)=Sc.AIC(5);
    bic_o(j)=Sc.BIC(1); bic_e(j)=Sc.BIC(2); bic_s(j)=Sc.BIC(3); bic_c(j)=Sc.BIC(4); bic_g(j)=Sc.BIC(5);
    ksp_o(j)=Sc.KSp(1); ksp_e(j)=Sc.KSp(2); ksp_s(j)=Sc.KSp(3); ksp_c(j)=Sc.KSp(4); ksp_g(j)=Sc.KSp(5);
end

%%

% Plot trailing fits.
figure(1); clf;
subplot(311);
loglog(Ts(2:end),1./diff(Ts),'ok','DisplayName','Data Sample'); hold on;
loglog(t,Sf.PDFo,':', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Regular Fit');
loglog(t,Sc.PDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Piecewise Fit');
%loglog(t,Sf.PDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
%loglog(t,Sf.PDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
%loglog(t,Sf.PDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
%loglog(t,Sf.PDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%loglog(bounds(1)*[1 1], ylim(),'--k');
%loglog(bounds(2)*[1 1], ylim(),'--k');
loglog(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
ylabel('Trailing Rate (events/day)'); xlabel('Time (days)');
legend('Location','southwest');
xlim(bounds);
subplot(312);
histogram(Ts,round(2*sqrt(length(Ts))),'Normalization','countdensity','DisplayName','Binned Data Sample'); hold on;
Ylim=ylim();
plot(t,Sf.PDFo,':', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Regular Fit');
plot(t,Sc.PDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Piecewise Fit');
%plot(t,Sf.PDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
%plot(t,Sf.PDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
%plot(t,Sf.PDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
%plot(t,Sf.PDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%plot(bounds(1)*[1 1], ylim(),'--k');
%plot(bounds(2)*[1 1], ylim(),'--k');
plot(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
ylabel('Trailing Rate (events/day)'); xlabel('Time (days)');
legend('Location','northeast');
xlim(bounds); ylim(Ylim*1.1);
subplot(313);
plot(Ts,1:length(Ts),'-k', 'LineWidth',2,'DisplayName','Data Sample'); hold on;
plot(t,Sf.CDFo,':', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Regular Fit');
plot(t,Sc.CDFo,'-', 'LineWidth',lw,'Color','#0000FF','DisplayName','Omori Piecewise Fit');
%plot(t,Sf.CDFe,'-', 'LineWidth',lw,'Color','#FF0000','DisplayName','Exponential Fit');
%plot(t,Sf.CDFs,'-', 'LineWidth',lw,'Color','#EDB120','DisplayName','Stretched Fit');
%plot(t,Sf.CDFc,'-', 'LineWidth',lw,'Color','#FF00FF','DisplayName','Cut-off Fit');
%plot(t,Sf.CDFg,'-', 'LineWidth',lw,'Color','#00FF00','DisplayName','Gamma Fit');
%plot(bounds(1)*[1 1], ylim(),'--k');
%plot(bounds(2)*[1 1], ylim(),'--k');
plot(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
ylabel('Trailing Counts'); xlabel('Time (days)');
legend('Location','southeast');
ylim([0 1.1]*length(Ts));
xlim(bounds);

% Plot performance metrics as a function of split-time.
figure(2); clf;
subplot(311);
plot(Tc,ll_o,'-o','Color','#0000FF','DisplayName','Omori'); hold on;
plot(Tc,ll_e,'-o','Color','#FF0000','DisplayName','Exponential');
plot(Tc,ll_s,'-o','Color','#EDB120','DisplayName','Stretched');
plot(Tc,ll_c,'-o','Color','#FF00FF','DisplayName','Cut-off');
plot(Tc,ll_g,'-o','Color','#00FF00','DisplayName','Gamma');
%plot(xlim(),Sf.LL(1)*[1 1],'-');
xlabel('Cut-off time (days)'); ylabel('Log-likelihood');
legend();
subplot(312);
plot(Tc,aic_o,'-o','Color','#0000FF','DisplayName','Omori'); hold on;
plot(Tc,aic_e,'-o','Color','#FF0000','DisplayName','Exponential');
plot(Tc,aic_s,'-o','Color','#EDB120','DisplayName','Stretched');
plot(Tc,aic_c,'-o','Color','#FF00FF','DisplayName','Cut-off');
plot(Tc,aic_g,'-o','Color','#00FF00','DisplayName','Gamma');
%plot(xlim(),Sf.AIC(1)*[1 1],'-');
xlabel('Cut-off time (days)'); ylabel('AIC');
subplot(313);
plot(Tc,bic_o,'-o','Color','#0000FF','DisplayName','Omori'); hold on;
plot(Tc,bic_e,'-o','Color','#FF0000','DisplayName','Exponential');
plot(Tc,bic_s,'-o','Color','#EDB120','DisplayName','Stretched');
plot(Tc,bic_c,'-o','Color','#FF00FF','DisplayName','Cut-off');
plot(Tc,bic_g,'-o','Color','#00FF00','DisplayName','Gamma');
%plot(xlim(),Sf.BIC(1)*[1 1],'-');
xlabel('Cut-off time (days)'); ylabel('BIC');

% Plot the CDF residuals.
figure(3); clf;
semilogx(t,zeros(size(t)),'-k'); hold on;
semilogx(t,residuals1(Ts,t,Sf.CDFo),':','Color','#0000FF','DisplayName','Omori Regular Fit');
semilogx(t,residuals1(Ts,t,Sc.CDFo),'-','Color','#0000FF','DisplayName','Omori Piecewise Fit');
plot(Tc(j)*[1 1], ylim(),':k','DisplayName','Optimal Changepoint');
plot(0.2648*1*[1 1],ylim(),'--k','DisplayName','Flowback Start'); % Start of flow back [digitized from Häring 2008].
plot(0.8553*1*[1 1],ylim(),'--k','DisplayName','Flowback Rate Decrease'); % first big decrease in flow back rate.
xlabel('Time (days)'); ylabel('Cumulative Density Residual (-)');
xlim([min(t) max(t)]);
legend('Location','southwest');


max(ll_o)
min(aic_o)
min(bic_o)
max(ksp_o)

Sf.LL(1)
Sf.AIC(1)
Sf.BIC(1)
Sf.KSp(1)




% Subroutines.
function res = residuals1(xd,xf,yf)
  % Local subroutine that calculates the CDF residuals (on fitted x-axis).
  
  yd=interp1(xd+1e-10*(1:length(xd))',1:length(xd),xf,'nearest','extrap');
  res=(yf-yd)/length(xd);
  
end