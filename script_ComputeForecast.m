% Script to compute the point-wise adaptive ensemble forecast.
clear;

% Predefine some values.
Nf=1e6;
Ns=5e3;
i1=10;
Forecast_case='SSFS-2005-d';
ic=1;

% Load the in the data structure and keep just the SSFS results.
load('DataFits.mat','Sf');
Sg=Sf(contains({Sf.ID},Forecast_case));

% Get the required information for the pointwise fits of stimulation D.
t=Sg(ic).t;
Ts=Sg(ic).Ts;
bounds_g=[0.9*min(Ts) max(Ts)+0.01];
bounds_f=bounds_g;

% Make the initial guesses for parameters, based on stimulations A-C.
%Po=[mean(arrayfun(@(S) S.Po(1), Sg(1:3))) mean(arrayfun(@(S) S.Po(2), Sg(1:3))) mean(arrayfun(@(S) S.Po(3), Sg(1:3)))];
%Pe=[mean(arrayfun(@(S) S.Pe(1), Sg(1:3))) mean(arrayfun(@(S) S.Pe(2), Sg(1:3)))];
%Ps=[mean(arrayfun(@(S) S.Ps(1), Sg(1:3))) mean(arrayfun(@(S) S.Ps(2), Sg(1:3))) mean(arrayfun(@(S) S.Ps(3), Sg(1:3))) mean(arrayfun(@(S) S.Ps(4), Sg(1:3)))];
%Pc=[mean(arrayfun(@(S) S.Pc(1), Sg(1:3))) mean(arrayfun(@(S) S.Pc(2), Sg(1:3))) mean(arrayfun(@(S) S.Pc(3), Sg(1:3))) mean(arrayfun(@(S) S.Pc(4), Sg(1:3)))];
%Pg=[mean(arrayfun(@(S) S.Pg(1), Sg(1:3))) mean(arrayfun(@(S) S.Pg(2), Sg(1:3))) mean(arrayfun(@(S) S.Pg(3), Sg(1:3))) mean(arrayfun(@(S) S.Pg(4), Sg(1:3)))];
Po=Sf(ic).Po; Pe=Sf(ic).Pe; Ps=Sf(ic).Ps; Pc=Sf(ic).Pc; Pg=Sf(ic).Pg;
%Po=[]; Pe=[]; Ps=[]; Pc=[]; Pg=[];

% Loop over each of the parts.
for i=i1:length(Ts)
    
    % Output to the screen.
    i
    length(Ts)
    
    % Get bounds for this iteration.
    bounds_f(2)=max(Ts(1:i))+0.01;
    
    % Fit.
    [sf]=Fit_Trailing(Ts(1:i),bounds_f, Nf,t, Po,Pe,Ps,Pc,Pg);
    for j=1:5
        [sf]=Fit_Trailing(Ts(1:i),bounds_f, Nf,t, sf.Po,sf.Pe,sf.Ps,sf.Pc,sf.Pg);
    end
    sf.i=i;
    
    % Update the initial guesses.
    Po=sf.Po;
    Pe=sf.Pe;
    Ps=sf.Ps;
    Pc=sf.Pc;
    Pg=sf.Pg;
    
    % Get some other metrics.
    % The expected number of events at end time.
    % The fit to the forecast part of the data.

    % Stuff the data into the structure.
    if(i==i1)
        Sf=sf;
    else
        Sf(end+1)=sf;
    end

end

% Save the data file.
save('DataForecast.mat','Sf','Sg','-v7.3');




% % Plot trailing fits.
% figure(2); clf;
% subplot(311);
% loglog(Ts(2:end),1./diff(Ts),'ok','DisplayName','Data Sample'); hold on;
% loglog(t,Sf(end).PDFo,'-b','DisplayName','Omori Fit');
% loglog(t,Sf(end).PDFe,'-r','DisplayName','Exponential Fit');
% loglog(t,Sf(end).PDFs,'-y','DisplayName','Stretched Fit');
% loglog(t,Sf(end).PDFc,'-m','DisplayName','Cut-off Fit');
% loglog(t,Sf(end).PDFg,'-g','DisplayName','Gamma Fit');
% loglog(bounds_f(1)*[1 1], ylim(),'--k');
% loglog(bounds_f(2)*[1 1], ylim(),'--k');
% ylabel('Aftershock Rate (events/day)'); xlabel('Time (days)');
% legend('Location','southwest');
% xlim(bounds_g);
% subplot(312);
% histogram(Ts,round(2*sqrt(length(Ts))),'Normalization','countdensity','DisplayName','Binned Data Sample'); hold on;
% Ylim=ylim();
% plot(t,Sf(end).PDFo,'-b','DisplayName','Omori Fit');
% plot(t,Sf(end).PDFe,'-r','DisplayName','Exponential Fit');
% plot(t,Sf(end).PDFs,'-y','DisplayName','Stretched Fit');
% plot(t,Sf(end).PDFc,'-m','DisplayName','Cut-off Fit');
% plot(t,Sf(end).PDFg,'-g','DisplayName','Gamma Fit');
% plot(bounds_f(1)*[1 1], ylim(),'--k');
% plot(bounds_f(2)*[1 1], ylim(),'--k');
% ylabel('Aftershock Rate (events/day)'); xlabel('Time (days)');
% legend('Location','northeast');
% xlim(bounds_g); ylim(Ylim*1.1);
% subplot(313);
% plot(Ts,1:length(Ts),'-k','DisplayName','Data Sample'); hold on;
% plot(t,Sf(end).CDFo,'-b','DisplayName','Omori Fit');
% plot(t,Sf(end).CDFe,'-r','DisplayName','Exponential Fit');
% plot(t,Sf(end).CDFs,'-y','DisplayName','Stretched Fit');
% plot(t,Sf(end).CDFc,'-m','DisplayName','Cut-off Fit');
% plot(t,Sf(end).CDFg,'-g','DisplayName','Gamma Fit');
% plot(bounds_f(1)*[1 1], ylim(),'--k');
% plot(bounds_f(2)*[1 1], ylim(),'--k');
% ylabel('Aftershock Counts'); xlabel('Time (days)');
% legend('Location','southeast');
% ylim([0 1.1]*length(Ts));
% xlim(bounds_g);
