% Script to make Figure 2.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Define some variables.
n=length(Sf);

% Factors for TDT, Bath's law.
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
Fs=Fs(Fs<50);
Fc=Fc(Fc<50);
Fg=Fg(Fg<50);

% Fit each of the factor distributions.
pHatO1=lognfit(Fo); pHatO2=expfit(Fo);
pHatE1=lognfit(Fe); pHatE2=expfit(Fe);
pHatS1=lognfit(Fs); pHatS2=expfit(Fs);
pHatC1=lognfit(Fc); pHatC2=expfit(Fc);
pHatG1=lognfit(Fg); pHatG2=expfit(Fg);
pHatA1=lognfit(F);  pHatA2=expfit(F);

% Plot Bath's law decay factors.
figure(2); clf;
histogram(F,round(5*sqrt(length(F))),'Normalization','pdf'); hold on;
plot(0:1:round(max(F)),lognpdf(0:1:round(max(F)),pHatA1(1),pHatA1(2)),'-b');
plot(0:1:round(max(F)),exppdf(0:1:round(max(F)),pHatA2),'-r');
xlabel('Median Typical Decay Time (days) [TDT: N^m(\infty)/n^m(t=0)] '); ylabel('Probability Density');
xlim([0 50]);

