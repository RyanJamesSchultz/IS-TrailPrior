% Script to fit and examine performance metrics & residuals for the synthetic-case data.
clear;

% Predefine some values.
dt=1e-3;
t=0:dt:1e3; t(1)=[];
Nsamp=1e4;
Nf=Nsamp*1e2;
n=1e1;

% Preallocate.
St=struct();

% Loop.
for i=1:n
    
    %%% Omori.
    K=1; c=0.01; p=1.2;
    initial=[Nsamp*(p-1)/(c^(1-p)),c,p];
    type_flag='Omori';
    guess_o=[1,0.01,1.2];
    guess_e=[1,3];
    guess_s=[1,0.01,0.005,0.15];
    guess_c=[1,0.01,450,0.90];
    guess_g=[1,100,1e-4,1.2];
    
    %%% Exponential.
    %n0=1; tau=2;
    %initial=[Nsamp/tau,tau];
    %type_flag='Exponential';
    %guess_o=[1,1,2];
    %guess_e=[1,3];
    %guess_s=[1,2,2,1];
    %guess_c=[1,0.002,2,0.001];
    %guess_g=[1,50,1e-6,1.1];
    
    %%% Stretched exponential.
    %Ns=1; t0=150; d=0.01; q=0.1;
    %initial=[Nsamp,t0,d,q];
    %type_flag='Stretched';
    %guess_o=[1,0.01,1.2];
    %guess_e=[1,100];
    %guess_s=[1,150,0.01,0.1];
    %guess_c=[1,1e-7,90,0.01];
    %guess_g=[1,50,1e-6,1.1];
    
    %%% Cut-off power law.
    %K=1; c=0.01; T=500; p=0.9; 
    %initial=[Nsamp/(T^(1-p)*exp(c/T)*gammainc(c/T,1-p,'upper')*gamma(1-p)),c,T,p];
    %type_flag='Cut-off';
    %guess_o=[1,0.01,1.1];
    %guess_e=[1,70];
    %guess_s=[1,20,0.005,0.25];
    %guess_c=[1,1e-7,90,0.01];
    %guess_g=[1,30,1e-6,1.1];
    
    %%% Gamma.
    %A=1; lb=50; la=1e-5; q=1.1;
    %initial=[Nsamp/((1/(q-1))*(lb^(q-1)-la^(q-1))),lb,la,q];
    %type_flag='Gamma';
    %guess_o=[1,0.01,1.2];
    %guess_e=[1,3];
    %guess_s=[1,0.01,0.01,0.1];
    %guess_c=[1,0.005,900,1.00];
    %guess_g=[1,50,1e-5,1.2];
    
    % Generate the test sample.
    [PDFt,CDFt,Nt]=EQ_Rate_Decay(t,type_flag,initial);
    Ts=EQ_Rate_Decay_Rand(t,[Nsamp 1],type_flag,initial);
    St(1,1,i).PDFt=PDFt; St(1,1,i).CDFt=CDFt; St(1,1,i).Ts=Ts;
    bounds=[0.9*min(Ts),1.1*max(Ts)];
    
    % Fit each model to the test sample.
    [s]=Fit_Trailing(Ts,bounds, Nf,t, guess_o,guess_e,guess_s,guess_c,guess_g);

    % Save the data to the structure.
    if(i==1)
        Sf=s;
    else
        Sf(end+1)=s;
    end
    
end


%%


% Report some values.
fprintf(1,'\n');
fprintf(1,'%s\n',['True ',type_flag]);
fprintf(1,'dAICc values\n');
fprintf(1,'%0.1e %0.1e %0.1e %0.1e %0.1e\n',mean(arrayfun(@(S) S.dAIC(1), Sf)),mean(arrayfun(@(S) S.dAIC(2), Sf)),mean(arrayfun(@(S) S.dAIC(3), Sf)),mean(arrayfun(@(S) S.dAIC(4), Sf)),mean(arrayfun(@(S) S.dAIC(5), Sf)));
fprintf(1,'dBIC values\n');
fprintf(1,'%0.1e %0.1e %0.1e %0.1e %0.1e\n',mean(arrayfun(@(S) S.dBIC(1), Sf)),mean(arrayfun(@(S) S.dBIC(2), Sf)),mean(arrayfun(@(S) S.dBIC(3), Sf)),mean(arrayfun(@(S) S.dBIC(4), Sf)),mean(arrayfun(@(S) S.dBIC(5), Sf)));
fprintf(1,'Unbinned R2 values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',mean(arrayfun(@(S) S.R2r(1), Sf)),mean(arrayfun(@(S) S.R2r(2), Sf)),mean(arrayfun(@(S) S.R2r(3), Sf)),mean(arrayfun(@(S) S.R2r(4), Sf)),mean(arrayfun(@(S) S.R2r(5), Sf)));
fprintf(1,'Binned R2 values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',mean(arrayfun(@(S) S.R2b(1), Sf)),mean(arrayfun(@(S) S.R2b(2), Sf)),mean(arrayfun(@(S) S.R2b(3), Sf)),mean(arrayfun(@(S) S.R2b(4), Sf)),mean(arrayfun(@(S) S.R2b(5), Sf)));
fprintf(1,'log_{10} KS values\n');
fprintf(1,'%+7.4f %+7.4f %+7.4f %+7.4f %+7.4f\n',mean(arrayfun(@(S) S.KSp(1), Sf)),mean(arrayfun(@(S) S.KSp(2), Sf)),mean(arrayfun(@(S) S.KSp(3), Sf)),mean(arrayfun(@(S) S.KSp(4), Sf)),mean(arrayfun(@(S) S.KSp(5), Sf)));
fprintf(1,'\nOrdering\n');
fprintf(1,'Omori, Exponential, Stretch, Cut-off, Gamma\n');

% Plot some fit metrics.
figure(1); clf;
%subplot(231);
%histogram(-arrayfun(@(S) S.LL(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
%histogram(-arrayfun(@(S) S.LL(2), Sf),round(sqrt(n)),'DisplayName','Exp');
%histogram(-arrayfun(@(S) S.LL(3), Sf),round(sqrt(n)),'DisplayName','Str');
%histogram(-arrayfun(@(S) S.LL(4), Sf),round(sqrt(n)),'DisplayName','Cut');
%histogram(-arrayfun(@(S) S.LL(5), Sf),round(sqrt(n)),'DisplayName','Gam');
%xlabel('-LL'); ylabel('Count');
%legend('Location','northeast');
%title(['True ',type_flag]);
subplot(221);
histogram(arrayfun(@(S) S.dAIC(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.dAIC(2), Sf),round(sqrt(n)),'DisplayName','Exp');
histogram(arrayfun(@(S) S.dAIC(3), Sf),round(sqrt(n)),'DisplayName','Str');
histogram(arrayfun(@(S) S.dAIC(4), Sf),round(sqrt(n)),'DisplayName','Cut');
histogram(arrayfun(@(S) S.dAIC(5), Sf),round(sqrt(n)),'DisplayName','Gam');
xlabel('\DeltaAIC'); ylabel('Count');
legend('Location','northeast');
%title(['True ',type_flag]);
subplot(222);
histogram(arrayfun(@(S) S.dBIC(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.dBIC(2), Sf),round(sqrt(n)),'DisplayName','Exp');
histogram(arrayfun(@(S) S.dBIC(3), Sf),round(sqrt(n)),'DisplayName','Str');
histogram(arrayfun(@(S) S.dBIC(4), Sf),round(sqrt(n)),'DisplayName','Cut');
histogram(arrayfun(@(S) S.dBIC(5), Sf),round(sqrt(n)),'DisplayName','Gam');
xlabel('\DeltaBIC'); ylabel('Count');
legend('Location','northeast');
%title(['True ',type_flag]);
subplot(223);
histogram(arrayfun(@(S) S.R2b(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.R2b(2), Sf),round(sqrt(n)),'DisplayName','Exp');
histogram(arrayfun(@(S) S.R2b(3), Sf),round(sqrt(n)),'DisplayName','Str');
histogram(arrayfun(@(S) S.R2b(4), Sf),round(sqrt(n)),'DisplayName','Cut');
histogram(arrayfun(@(S) S.R2b(5), Sf),round(sqrt(n)),'DisplayName','Gam');
xlabel('Goodness-of-fit R^2 (binned)'); ylabel('Count');
legend('Location','northwest');
%title(['True ',type_flag]);
%subplot(235);
%histogram(arrayfun(@(S) S.R2r(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
%histogram(arrayfun(@(S) S.R2r(2), Sf),round(sqrt(n)),'DisplayName','Exp');
%histogram(arrayfun(@(S) S.R2r(3), Sf),round(sqrt(n)),'DisplayName','Str');
%histogram(arrayfun(@(S) S.R2r(4), Sf),round(sqrt(n)),'DisplayName','Cut');
%histogram(arrayfun(@(S) S.R2r(5), Sf),round(sqrt(n)),'DisplayName','Gam');
%xlabel('Goodness-of-fit R^2 (unbinned)'); ylabel('Count');
%legend('Location','northwest');
%title(['True ',type_flag]);
subplot(224);
histogram(arrayfun(@(S) S.KSp(1), Sf),round(sqrt(n)),'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.KSp(2), Sf),round(sqrt(n)),'DisplayName','Exp');
histogram(arrayfun(@(S) S.KSp(3), Sf),round(sqrt(n)),'DisplayName','Str');
histogram(arrayfun(@(S) S.KSp(4), Sf),round(sqrt(n)),'DisplayName','Cut');
histogram(arrayfun(@(S) S.KSp(5), Sf),round(sqrt(n)),'DisplayName','Gam');
xlabel('log_{10}KSp'); ylabel('Count');
legend('Location','northwest');
%title(['True ',type_flag]);

% Plot the data best fits.
figure(2); clf;
subplot(221);
loglog(St(1,1,i).Ts(2:end),1./diff(St(1,1,i).Ts),'ok','DisplayName','Synthetically Sampled'); hold on;
for i=1:n
    loglog(t,Sf(i).PDFo,'-', 'Color','#0000FF','DisplayName','Omori Fit');
    loglog(t,Sf(i).PDFe,'-', 'Color','#FF0000','DisplayName','Exponential Fit');
    loglog(t,Sf(i).PDFs,'-', 'Color','#EDB120','DisplayName','Stretched Fit');
    loglog(t,Sf(i).PDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    loglog(t,Sf(i).PDFg,'-', 'Color','#77AC30','DisplayName','Gamma Fit');
end
loglog(t,St(1,1,i).PDFt,'--k','LineWidth',2, 'DisplayName','Omori Real Model');
xlabel('Time (days)'); ylabel('Aftershock Rate (events/day)');
legend({'Synthetically Sampled','Omori Fit','Exponential Fit','Stretched Fit','Cut-off Fit','Gamma Fit'},'Location','southwest');
ylim([0.5*min(St(1,1,i).PDFt) 1.5*max(St(1,1,i).PDFt)]);
xlim([min(t) max(t)]);

subplot(222);
semilogx(t,zeros(size(t)),'-k', 'DisplayName','Zero Residual'); hold on;
for i=1:n
    semilogx(t,(Sf(i).PDFo-St(1,1,i).PDFt),'-', 'Color','#0000FF','DisplayName','Omori Fit');
    semilogx(t,(Sf(i).PDFe-St(1,1,i).PDFt),'-', 'Color','#FF0000','DisplayName','Exponential Fit');
    semilogx(t,(Sf(i).PDFs-St(1,1,i).PDFt),'-', 'Color','#EDB120','DisplayName','Stretched Fit');
    semilogx(t,(Sf(i).PDFc-St(1,1,i).PDFt),'-', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    semilogx(t,(Sf(i).PDFg-St(1,1,i).PDFt),'-', 'Color','#77AC30','DisplayName','Gamma Fit');
end
xlabel('Time (days)'); ylabel('Aftershock Rate Residuals (events/day)');
legend({'Zero Residual','Omori Fit','Exponential Fit','Stretched Fit','Cut-off Fit','Gamma Fit'},'Location','southeast');
ylim(2*Nsamp*[-1 +1]);
xlim([min(t) max(t)]);

subplot(223);
for i=1:n
    semilogx(sort(St(1,1,i).Ts),1:length(St(1,1,i).Ts),'-k','DisplayName','Synthetically Sampled'); hold on;
    semilogx(t,Sf(i).CDFo,'-b', 'Color','#0000FF','DisplayName','Omori Fit');
    semilogx(t,Sf(i).CDFe,'-r', 'Color','#FF0000','DisplayName','Exponential Fit');
    semilogx(t,Sf(i).CDFs,'-y', 'Color','#EDB120','DisplayName','Stretched Fit');
    semilogx(t,Sf(i).CDFc,'-m', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    semilogx(t,Sf(i).CDFg,'-g', 'Color','#77AC30','DisplayName','Gamma Fit');
end
ylabel('Aftershock Counts'); xlabel('Time (days)');
ylim([0 1.1]*length(St(1,1,i).Ts));
legend({'Synthetically Sampled','Omori Fit','Exponential Fit','Stretched Fit','Cut-off Fit','Gamma Fit'},'Location','northwest');
xlim([min(t) max(t)]);

subplot(224);
semilogx(t,zeros(size(t)),'-k','LineWidth',2,'DisplayName','Zero Residuals'); hold on;
for i=1:n
    semilogx(St(1,1,i).Ts,interp1(t,Sf(i).CDFo,St(1,1,i).Ts,'linear')-(1:length(St(1,1,i).Ts))','-b', 'Color','#0000FF','DisplayName','Omori Fit');
    semilogx(St(1,1,i).Ts,interp1(t,Sf(i).CDFe,St(1,1,i).Ts,'linear')-(1:length(St(1,1,i).Ts))','-r', 'Color','#FF0000','DisplayName','Exponential Fit');
    semilogx(St(1,1,i).Ts,interp1(t,Sf(i).CDFs,St(1,1,i).Ts,'linear')-(1:length(St(1,1,i).Ts))','-y', 'Color','#EDB120','DisplayName','Stretched Fit');
    semilogx(St(1,1,i).Ts,interp1(t,Sf(i).CDFc,St(1,1,i).Ts,'linear')-(1:length(St(1,1,i).Ts))','-m', 'Color','#FF00FF','DisplayName','Cut-off Fit');
    semilogx(St(1,1,i).Ts,interp1(t,Sf(i).CDFg,St(1,1,i).Ts,'linear')-(1:length(St(1,1,i).Ts))','-g', 'Color','#77AC30','DisplayName','Gamma Fit');
end
ylabel('Aftershock Count Residuals'); xlabel('Time (days)');
legend({'Zero Residuals','Omori Fit','Exponential Fit','Stretched Fit','Cut-off Fit','Gamma Fit'},'Location','northwest');
ylim(0.2*Nsamp*[-1 +1]);
xlim([min(t) max(t)]);

