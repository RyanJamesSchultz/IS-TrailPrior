% Script to make Figure S10.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Factors for Bath's law.
Fo=zeros(size(Sf)); Fe=Fo; Fs=Fo; Fc=Fo; Fg=Fo;
for i=1:length(Sf)
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Omori',Sf(i).Po);       Fo(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Exponential',Sf(i).Pe); Fe(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Stretched',Sf(i).Ps);   Fs(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Cut-off',Sf(i).Pc);     Fc(i)=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Gamma',Sf(i).Pg);       Fg(i)=Nt/nt;
end
F=median([Fo;Fe;Fs;Fc;Fg]);
F=F(F<100);
Fo=Fo(Fo<100);
Fe=Fe(Fe<100);
Fs=Fs(Fs<100);
Fc=Fc(Fc<100);
Fg=Fg(Fg<100);

% Fit each of the factor distributions.
pHatO1=lognfit(Fo); pHatO2=expfit(Fo);
pHatE1=lognfit(Fe); pHatE2=expfit(Fe);
pHatS1=lognfit(Fs); pHatS2=expfit(Fs);
pHatC1=lognfit(Fc); pHatC2=expfit(Fc);
pHatG1=lognfit(Fg); pHatG2=expfit(Fg);
pHatA1=lognfit(F);  pHatA2=expfit(F);

% Plot Bath's law decay factors.
figure(510); clf;
subplot(611);
histogram(Fo,round(5*sqrt(length(Fo))),'Normalization','pdf'); hold on;
plot(0:1:(max(Fo)),lognpdf(0:1:(max(Fo)),pHatO1(1),pHatO1(2)),'-b');
plot(0:1:(max(Fo)),exppdf(0:1:(max(Fo)),pHatO2),'-r');
xlabel('Omori Typical Decay Time (days) [TDT: N^o(\infty)/n^o(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);
subplot(612);
histogram(Fe,round(5*sqrt(length(Fe))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fe)),lognpdf(0:1:round(max(Fe)),pHatE1(1),pHatE1(2)),'-b');
plot(0:1:round(max(Fe)),exppdf(0:1:round(max(Fe)),pHatE2),'-r');
xlabel('Exponential Typical Decay Time (days) [TDT: N^e(\infty)/n^e(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);
subplot(613);
histogram(Fs,round(5*sqrt(length(Fs))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fs)),lognpdf(0:1:round(max(Fs)),pHatS1(1),pHatS1(2)),'-b');
plot(0:1:round(max(Fs)),exppdf(0:1:round(max(Fs)),pHatS2),'-r');
xlabel('Stretch Typical Decay Time (days) [TDT: N^s(\infty)/n^s(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);
subplot(614);
histogram(Fc,round(5*sqrt(length(Fc))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fc)),lognpdf(0:1:round(max(Fc)),pHatC1(1),pHatC1(2)),'-b');
plot(0:1:round(max(Fc)),exppdf(0:1:round(max(Fc)),pHatC2),'-r');
xlabel('Cut-off Typical Decay Time (days) [TDT: N^c(\infty)/n^c(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);
subplot(615);
histogram(Fg,round(5*sqrt(length(Fg))),'Normalization','pdf'); hold on;
plot(0:1:round(max(Fg)),lognpdf(0:1:round(max(Fg)),pHatG1(1),pHatG1(2)),'-b');
plot(0:1:round(max(Fg)),exppdf(0:1:round(max(Fg)),pHatG2),'-r');
xlabel('Gamma Typical Decay Time (days) [TDT: N^\gamma(\infty)/n^\gamma(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);
subplot(616);
histogram(F,round(5*sqrt(length(F))),'Normalization','pdf'); hold on;
plot(0:1:round(max(F)),lognpdf(0:1:round(max(F)),pHatA1(1),pHatA1(2)),'-b');
plot(0:1:round(max(F)),exppdf(0:1:round(max(F)),pHatA2),'-r');
xlabel('Median Typical Decay Time (days) [TDT: N^m(\infty)/n^m(t=0)] '); ylabel('Probability Density');
xlim([0 50]); ylim([0 0.6]);

