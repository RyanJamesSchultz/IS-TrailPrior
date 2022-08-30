% Script to examine the performance metrics for the fitted real-case data.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Define some variables.
n=length(Sf);

% Weights from AIC & BIC metrics.
Wa_av=[0 0 0 0 0];
Wb_av=Wa_av;
Wa_50=Wa_av;
Wb_50=Wa_av;
R2_av=Wa_av;
R2_50=Wa_av;
KS_av=Wa_av;
KS_50=Wa_av;
for i=1:length(Wa_av)
    Wa_av(i)=mean(arrayfun(@(S) S.Waic(i), Sf));
    Wb_av(i)=mean(arrayfun(@(S) S.Wbic(i), Sf));
    Wa_50(i)=median(arrayfun(@(S) S.Waic(i), Sf));
    Wb_50(i)=median(arrayfun(@(S) S.Wbic(i), Sf));
    R2_av(i)=mean(arrayfun(@(S) S.R2b(i), Sf));
    R2_50(i)=median(arrayfun(@(S) S.R2b(i), Sf));
    KS_av(i)=mean(arrayfun(@(S) S.KSp(i), Sf));
    KS_50(i)=median(arrayfun(@(S) S.KSp(i), Sf));
end
Wa_av=Wa_av/sum(Wa_av);
Wb_av=Wb_av/sum(Wb_av);
Wa_50=Wa_50/sum(Wa_50);
Wb_50=Wb_50/sum(Wb_50);

% Aggregate weights for logic trees.
W_av=mean([Wa_av;Wb_av]);
W_av=W_av/sum(W_av);
W_50=mean([Wa_50;Wb_50]);
W_50=W_50/sum(W_50);
W=mean([W_av;W_50]);

% TDT factors for Bath's law.
Fo=zeros(size(Sf)); Fe=Fo; Fs=Fo; Fc=Fo; Fg=Fo;
for i=1:length(Sf)
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Omori',Sf(i).Po);       Fo(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Exponential',Sf(i).Pe); Fe(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Stretched',Sf(i).Ps);   Fs(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Cut-off',Sf(i).Pc);     Fc(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Gamma',Sf(i).Pg);       Fg(i)=Nt/nt;
end
F=median([Fo;Fe;Fs;Fc;Fg]);
F(55)=Fs(55);
%F=F(F<1000);
Fo=Fo(Fo<100);
Fe=Fe(Fe<100);
Fs=Fs(Fs<50);
Fc=Fc(Fc<50);
Fg=Fg(Fg<50);

% Best KS scores.
KSbest=zeros(size(Sf));
for i =1:length(Sf)
    KSbest(i)=max(Sf(i).KSp);
end
length(find(KSbest>log10(0.05)))/length(KSbest)

% Fit each of the factor distributions.
pHatO1=lognfit(Fo); pHatO2=expfit(Fo);
pHatE1=lognfit(Fe); pHatE2=expfit(Fe);
pHatS1=lognfit(Fs); pHatS2=expfit(Fs);
pHatC1=lognfit(Fc); pHatC2=expfit(Fc);
pHatG1=lognfit(Fg); pHatG2=expfit(Fg);
pHatA1=lognfit(F);  pHatA2=expfit(F);



% Plot some fit metrics.
figure(1); clf;
subplot(231);
histogram(arrayfun(@(S) (S.Waic(1)+S.Wbic(1))/2, Sf),0:0.01:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) (S.Waic(2)+S.Wbic(2))/2, Sf),0:0.01:1,'DisplayName','Exp');
histogram(arrayfun(@(S) (S.Waic(3)+S.Wbic(3))/2, Sf),0:0.01:1,'DisplayName','Str');
histogram(arrayfun(@(S) (S.Waic(4)+S.Wbic(4))/2, Sf),0:0.01:1,'DisplayName','Cut');
histogram(arrayfun(@(S) (S.Waic(5)+S.Wbic(5))/2, Sf),0:0.01:1,'DisplayName','Gam');
xlabel('Wavg'); ylabel('Count');
legend('Location','northeast');
subplot(232);
histogram(arrayfun(@(S) S.Waic(1), Sf),0:0.01:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.Waic(2), Sf),0:0.01:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.Waic(3), Sf),0:0.01:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.Waic(4), Sf),0:0.01:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.Waic(5), Sf),0:0.01:1,'DisplayName','Gam');
xlabel('AIC Weight'); ylabel('Count');
legend('Location','northeast');
subplot(233);
histogram(arrayfun(@(S) S.Wbic(1), Sf),0:0.01:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.Wbic(2), Sf),0:0.01:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.Wbic(3), Sf),0:0.01:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.Wbic(4), Sf),0:0.01:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.Wbic(5), Sf),0:0.01:1,'DisplayName','Gam');
xlabel('BIC Weight'); ylabel('Count');
legend('Location','northeast');
subplot(234);
histogram(arrayfun(@(S) S.R2b(1), Sf),-4:0.01:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.R2b(2), Sf),-4:0.01:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.R2b(3), Sf),-4:0.01:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.R2b(4), Sf),-4:0.01:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.R2b(5), Sf),-4:0.01:1,'DisplayName','Gam');
xlabel('Goodness-of-fit R^2 (binned)'); ylabel('Count');
legend('Location','northwest');
subplot(235);
histogram(arrayfun(@(S) S.R2r(1), Sf),-4:0.01:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.R2r(2), Sf),-4:0.01:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.R2r(3), Sf),-4:0.01:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.R2r(4), Sf),-4:0.01:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.R2r(5), Sf),-4:0.01:1,'DisplayName','Gam');
xlabel('Goodness-of-fit R^2 (unbinned)'); ylabel('Count');
legend('Location','northwest');
subplot(236);
histogram(arrayfun(@(S) S.KSp(1), Sf),-10:0.1:0,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.KSp(2), Sf),-10:0.1:0,'DisplayName','Exp');
histogram(arrayfun(@(S) S.KSp(3), Sf),-10:0.1:0,'DisplayName','Str');
histogram(arrayfun(@(S) S.KSp(4), Sf),-10:0.1:0,'DisplayName','Cut');
histogram(arrayfun(@(S) S.KSp(5), Sf),-10:0.1:0,'DisplayName','Gam');
plot(log10(0.05)*[1 1], ylim(),'--k');
xlabel('logKS'); ylabel('Count');
legend('Location','northwest');

% Plot best KS-test p-values.
figure(3); clf;
histogram(KSbest,200); hold on;
plot(log10(0.05)*[1 1], ylim(),'--k');
xlabel('logKS'); ylabel('Count');

% Plot Bath's law decay factors.
figure(2); clf;
subplot(611);
histogram(Fo,round(5*sqrt(length(Fo))),'Normalization','pdf'); hold on;
plot(0:1:(max(Fo)),lognpdf(0:1:(max(Fo)),pHatO1(1),pHatO1(2)),'-b');
plot(0:1:(max(Fo)),exppdf(0:1:(max(Fo)),pHatO2),'-r');
xlabel('Omori Bath Factor [c/(p-1)] (days)'); ylabel('Count');
xlim([0 50]);
subplot(612);
histogram(Fe,round(5*sqrt(length(Fe))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fe)),lognpdf(0:1:round(max(Fe)),pHatE1(1),pHatE1(2)),'-b');
plot(0:1:round(max(Fe)),exppdf(0:1:round(max(Fe)),pHatE2),'-r');
xlabel('Exponential Mean Relataxation Time \tau (days)'); ylabel('Count');
xlim([0 50]);
subplot(613);
histogram(Fs,round(5*sqrt(length(Fs))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fs)),lognpdf(0:1:round(max(Fs)),pHatS1(1),pHatS1(2)),'-b');
plot(0:1:round(max(Fs)),exppdf(0:1:round(max(Fs)),pHatS2),'-r');
xlabel('Stretch Bath Factor (days)'); ylabel('Count');
xlim([0 50]);
subplot(614);
histogram(Fc,round(5*sqrt(length(Fc))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fc)),lognpdf(0:1:round(max(Fc)),pHatC1(1),pHatC1(2)),'-b');
plot(0:1:round(max(Fc)),exppdf(0:1:round(max(Fc)),pHatC2),'-r');
xlabel('Cut-off Bath Factor (days)'); ylabel('Count');
xlim([0 50]);
subplot(615);
histogram(Fg,round(5*sqrt(length(Fg))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fg)),lognpdf(0:1:round(max(Fg)),pHatG1(1),pHatG1(2)),'-b');
plot(0:1:round(max(Fg)),exppdf(0:1:round(max(Fg)),pHatG2),'-r');
xlabel('Cut-off Bath Factor (days)'); ylabel('Count');
xlim([0 50]);
subplot(616);
histogram(F,round(5*sqrt(length(F))),'Normalization','pdf'); hold on;
plot(0:1:round(max(F)),lognpdf(0:1:round(max(F)),pHatA1(1),pHatA1(2)),'-b');
plot(0:1:round(max(F)),exppdf(0:1:round(max(F)),pHatA2),'-r');
xlabel('Aggregate Bath Factor (days)'); ylabel('Count');
xlim([0 50]);

