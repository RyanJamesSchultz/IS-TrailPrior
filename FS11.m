% Script to make Figure S11.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Plot some fit metrics.
figure(511); clf;
subplot(221);
histogram(arrayfun(@(S) S.Waic(1), Sf),0:0.025:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.Waic(2), Sf),0:0.025:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.Waic(3), Sf),0:0.025:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.Waic(4), Sf),0:0.025:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.Waic(5), Sf),0:0.025:1,'DisplayName','Gam');
xlabel('W^{AIC}'); ylabel('Count');
legend('Location','northeast');
xlim([0 1]);
subplot(222);
histogram(arrayfun(@(S) S.Wbic(1), Sf),0:0.025:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.Wbic(2), Sf),0:0.025:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.Wbic(3), Sf),0:0.025:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.Wbic(4), Sf),0:0.025:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.Wbic(5), Sf),0:0.025:1,'DisplayName','Gam');
xlabel('W^{BIC}'); ylabel('Count');
legend('Location','northeast');
xlim([0 1]);
subplot(223);
histogram(arrayfun(@(S) S.R2b(1), Sf),-4:0.025:1,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.R2b(2), Sf),-4:0.025:1,'DisplayName','Exp');
histogram(arrayfun(@(S) S.R2b(3), Sf),-4:0.025:1,'DisplayName','Str');
histogram(arrayfun(@(S) S.R2b(4), Sf),-4:0.025:1,'DisplayName','Cut');
histogram(arrayfun(@(S) S.R2b(5), Sf),-4:0.025:1,'DisplayName','Gam');
xlabel('R^2'); ylabel('Count');
legend('Location','northwest');
xlim([0 1]);
subplot(224);
histogram(arrayfun(@(S) S.KSp(1), Sf),-10:0.25:0,'DisplayName','Omori'); hold on;
histogram(arrayfun(@(S) S.KSp(2), Sf),-10:0.25:0,'DisplayName','Exp');
histogram(arrayfun(@(S) S.KSp(3), Sf),-10:0.25:0,'DisplayName','Str');
histogram(arrayfun(@(S) S.KSp(4), Sf),-10:0.25:0,'DisplayName','Cut');
histogram(arrayfun(@(S) S.KSp(5), Sf),-10:0.25:0,'DisplayName','Gam');
plot(log10(0.05)*[1 1], ylim(),'--k','DisplayName','5% threshold');
xlabel('log_{10}KSp'); ylabel('Count');
legend('Location','northwest');
xlim([-10 0]);
