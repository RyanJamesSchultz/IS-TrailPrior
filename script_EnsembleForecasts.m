% Script to plot the ensemble forecasting results.  For analysis and QC'ing.
clear;

% Load the in the data structures for plotting.
load('DataForecastSSFS2005.mat','Sf','Sg');

% Get the required information for the pointwise fits of stimulation D.
t=Sg(end).t;
Ts=Sg(end).Ts;
bounds_g=[0.9*min(Ts) 1.1*max(Ts)];
bounds_f=bounds_g;

% Make the initial guesses for parameters, based on stimulations A-C.
Po_a=[mean(arrayfun(@(S) S.Po(1), Sg(end))) mean(arrayfun(@(S) S.Po(2), Sg(end))) mean(arrayfun(@(S) S.Po(3), Sg(end)))];
Pe_a=[mean(arrayfun(@(S) S.Pe(1), Sg(end))) mean(arrayfun(@(S) S.Pe(2), Sg(end)))];
Ps_a=[mean(arrayfun(@(S) S.Ps(1), Sg(end))) mean(arrayfun(@(S) S.Ps(2), Sg(end))) mean(arrayfun(@(S) S.Ps(3), Sg(end))) mean(arrayfun(@(S) S.Ps(4), Sg(end)))];
Pc_a=[mean(arrayfun(@(S) S.Pc(1), Sg(end))) mean(arrayfun(@(S) S.Pc(2), Sg(end))) mean(arrayfun(@(S) S.Pc(3), Sg(end))) mean(arrayfun(@(S) S.Pc(4), Sg(end)))];
Pg_a=[mean(arrayfun(@(S) S.Pg(1), Sg(end))) mean(arrayfun(@(S) S.Pg(2), Sg(end))) mean(arrayfun(@(S) S.Pg(3), Sg(end))) mean(arrayfun(@(S) S.Pg(4), Sg(end)))];
N_a=mean(arrayfun(@(S) length(S.Ts), Sg(end)));
N_a=10;

% Make the initial guesses for the weights.
W=mean([0.2751 0.1941 0.1974 0.1016 0.2318; 0.3382 0.2207 0.1556 0.0918 0.1937]); W=W/sum(W);

% Loop over each of the parts.
for i=1:length(Sf)
    
    % Get bounds for this iteration.
    bounds_f(2)=max(Ts(1:Sf(i).i))+1e-2;
    
    % Get new parameter guesses.
    N_o=Sf(i).i;
    Po_p=(N_o*Sf(i).Po+N_a*Po_a)/(N_a+N_o); Po_p(1)=Sf(i).Po(1);
    Pe_p=(N_o*Sf(i).Pe+N_a*Pe_a)/(N_a+N_o); Pe_p(1)=Sf(i).Pe(1);
    Ps_p=(N_o*Sf(i).Ps+N_a*Ps_a)/(N_a+N_o); Ps_p(1)=Sf(i).Ps(1);
    Pc_p=(N_o*Sf(i).Pc+N_a*Pc_a)/(N_a+N_o); Pc_p(1)=Sf(i).Pc(1);
    Pg_p=(N_o*Sf(i).Pg+N_a*Pg_a)/(N_a+N_o); Pg_p(1)=Sf(i).Pg(1);
    
    % Get the new predictions.
    [PDFo,CDFo,No]=EQ_Rate_Decay(t,'Omori',Po_p);
    [PDFe,CDFe,Ne]=EQ_Rate_Decay(t,'Exponential',Pe_p);
    [PDFs,CDFs,Ns]=EQ_Rate_Decay(t,'Stretched',Ps_p);
    [PDFc,CDFc,Nc]=EQ_Rate_Decay(t,'Cut-off',Pc_p);
    [PDFg,CDFg,Ng]=EQ_Rate_Decay(t,'Gamma',Pg_p);
    
    % Get the ensemble weights and prediction.
    W=mean([Sf(i).Wbic;Sf(i).Waic]); W=W/sum(W);
    PDFa=sum([W(1).*Sf(i).PDFo; W(2).*Sf(i).PDFe; W(3).*Sf(i).PDFs; W(4).*Sf(i).PDFc; W(5).*Sf(i).PDFg;]);
    CDFa=sum([W(1).*Sf(i).CDFo; W(2).*Sf(i).CDFe; W(3).*Sf(i).CDFs; W(4).*Sf(i).CDFc; W(5).*Sf(i).CDFg;]);
    
    % Get some other metrics.
    % The expected number of events at end time.
    % The fit to the forecast part of the data.
    
    
    
    % Plot trailing fits.
    figure(1); clf;
    subplot(311);
    loglog(Ts(2:end),1./diff(Ts),'ok','DisplayName','Data Sample'); hold on;
    loglog(t,Sf(i).PDFo,'-', 'Color','#0000FF','DisplayName','Omori Fit');
    loglog(t,Sf(i).PDFe,'-', 'Color','#FF0000','DisplayName','Exponential Fit');
    loglog(t,Sf(i).PDFs,'-', 'Color','#EDB120','DisplayName','Stretched Fit');
    loglog(t,Sf(i).PDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    loglog(t,Sf(i).PDFg,'-', 'Color','#77AC30','DisplayName','Gamma Fit');
    loglog(t,PDFa,'-c','DisplayName','Ensemble Forecast');
    %loglog(t,PDFo,'-b','DisplayName','Omori Posterior');
    %loglog(t,PDFe,'-r','DisplayName','Exponential Posterior');
    %loglog(t,PDFs,'-y','DisplayName','Stretched Posterior');
    %loglog(t,PDFc,'-m','DisplayName','Cut-off Posterior');
    %loglog(t,PDFg,'-g','DisplayName','Gamma Posterior');
    loglog(bounds_f(1)*[1 1], ylim(),'--k');
    loglog(bounds_f(2)*[1 1], ylim(),'--k');
    ylabel('Aftershock Rate (events/day)'); xlabel('Time (days)');
    legend('Location','southwest');
    xlim(bounds_g); ylim([0.9*min(1./diff(Ts)) 1.1*max(1./diff(Ts))]);
    subplot(312);
    histogram(Ts,round(2*sqrt(length(Ts))),'Normalization','countdensity','DisplayName','Binned Data Sample'); hold on;
    Ylim=ylim();
    plot(t,Sf(i).PDFo,'-', 'Color','#0000FF','DisplayName','Omori Fit');
    plot(t,Sf(i).PDFe,'-', 'Color','#FF0000','DisplayName','Exponential Fit');
    plot(t,Sf(i).PDFs,'-', 'Color','#EDB120','DisplayName','Stretched Fit');
    plot(t,Sf(i).PDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    plot(t,Sf(i).PDFg,'-', 'Color','#77AC30','DisplayName','Gamma Fit');
    plot(t,PDFa,'-c','DisplayName','Ensemble Forecast');
    %plot(t,PDFo,'-b','DisplayName','Omori Posterior');
    %plot(t,PDFe,'-r','DisplayName','Exponential Posterior');
    %plot(t,PDFs,'-y','DisplayName','Stretched Posterior');
    %plot(t,PDFc,'-m','DisplayName','Cut-off Posterior');
    %plot(t,PDFg,'-g','DisplayName','Gamma Posterior');
    plot(bounds_f(1)*[1 1], ylim(),'--k');
    plot(bounds_f(2)*[1 1], ylim(),'--k');
    ylabel('Aftershock Rate (events/day)'); xlabel('Time (days)');
    legend('Location','northeast');
    xlim(bounds_g); ylim(Ylim*1.1);
    subplot(313);
    plot(Ts,1:length(Ts),'-k','DisplayName','Data Sample'); hold on;
    plot(t,Sf(i).CDFo,'-', 'Color','#0000FF','DisplayName','Omori Fit');
    plot(t,Sf(i).CDFe,'-', 'Color','#FF0000','DisplayName','Exponential Fit');
    plot(t,Sf(i).CDFs,'-', 'Color','#EDB120','DisplayName','Stretched Fit');
    plot(t,Sf(i).CDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    plot(t,Sf(i).CDFg,'-', 'Color','#77AC30','DisplayName','Gamma Fit');
    plot(t,CDFa,'-c','DisplayName','Ensemble Forecast');
    %plot(t,CDFo,'-b','DisplayName','Omori Posterior');
    %plot(t,CDFe,'-r','DisplayName','Exponential Posterior');
    %plot(t,CDFs,'-y','DisplayName','Stretched Posterior');
    %plot(t,CDFc,'-m','DisplayName','Cut-off Posterior');
    %plot(t,CDFg,'-g','DisplayName','Gamma Posterior');
    plot(bounds_f(1)*[1 1], ylim(),'--k');
    plot(bounds_f(2)*[1 1], ylim(),'--k');
    plot(2*bounds_f(2)*[1 1], ylim(),':k');
    ylabel('Aftershock Counts'); xlabel('Time (days)');
    legend('Location','southeast');
    ylim([0 1.1]*length(Ts));
    xlim(bounds_g);
    
    i+Sf(1).i
    W
    pause;
    
end


