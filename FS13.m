% Script to make figure S13.
clear;

% Load in a predefined structure for convenience.
tableTData;

% Predefine some values.
GREY=[0.85,0.85,0.85];
Nf=1e6;
ID='Basel';
lw=1;

% Define the changepoint time.
Tc=2:0.5:15;
%Tc=(9:0.1:10)+0.01;
%Tc=[9.77:0.01:9.79,9.81:0.01:9.83];
%Tc=9.81
%Tc=3.5:0.1:5.0;
%Tc=4.0:0.01:4.2;
%Tc=4.08;

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






% Plot performance metrics as a function of split-time.
figure(513); clf;
subplot(211);
plot(Tc,aic_o,'-o','Color','#0000FF','DisplayName','Omori'); hold on;
plot(Tc,aic_e,'-o','Color','#FF0000','DisplayName','Exponential');
plot(Tc,aic_s,'-o','Color','#EDB120','DisplayName','Stretched');
plot(Tc,aic_c,'-o','Color','#FF00FF','DisplayName','Cut-off');
plot(Tc,aic_g,'-o','Color','#00FF00','DisplayName','Gamma');
xlabel('Split-time (days)'); ylabel('AIC');
xlim([min(Tc) max(Tc)]); ylim([-3210 -3140]);
legend();
subplot(212);
plot(Tc,bic_o,'-o','Color','#0000FF','DisplayName','Omori'); hold on;
plot(Tc,bic_e,'-o','Color','#FF0000','DisplayName','Exponential');
plot(Tc,bic_s,'-o','Color','#EDB120','DisplayName','Stretched');
plot(Tc,bic_c,'-o','Color','#FF00FF','DisplayName','Cut-off');
plot(Tc,bic_g,'-o','Color','#00FF00','DisplayName','Gamma');
xlabel('Split-time (days)'); ylabel('BIC');
xlim([min(Tc) max(Tc)]); ylim([-3180 -3120]);
legend();


