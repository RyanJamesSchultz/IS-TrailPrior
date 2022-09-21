% Script to make Figure S12.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Best KS scores.
KSbest=zeros(size(Sf));
for i =1:length(Sf)
    KSbest(i)=max(Sf(i).KSp);
end
length(find(KSbest>log10(0.05)))/length(KSbest)
length(find(KSbest<=log10(0.05)))

% Plot best KS-test p-values.
figure(3); clf;
histogram(KSbest,200); hold on;
plot(log10(0.05)*[1 1], ylim(),'--k');
xlabel('log_{10}KSp'); ylabel('Count');

