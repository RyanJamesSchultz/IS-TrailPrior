% Script to make Figure 4.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Truncate data to just the halted cases.
i=[17 19 20 21 46]; % 5% confidence.
%i=[8 14 16 17 19 20 21 22 26 46]; % 10% confidence.
Sf=Sf(i);

% Print out case names.
Sf.ID

% Plot the CDF residuals.
figure(1); clf;
subplot(121);
semilogx([1e-3 1e3],[0 0],'-k'); hold on;
% Loop over every case study.
for i=1:length(Sf)

    % Get the value of TDT.
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Omori',Sf(i).Po); Fto=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Exponential',Sf(i).Pe); Fte=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Stretched',Sf(i).Ps); Fts=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Cut-off',Sf(i).Pc); Ftc=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Gamma',Sf(i).Pg); Ftg=Nt/nt;
    Ft=median([Fto Fte Fts Ftc Ftg]);
    
    % Plot the residuals of each model.
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFo),'-', 'Color','#0000FF','DisplayName','Omori');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFe),'-', 'Color','#FF0000','DisplayName','Exp');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFs),'-', 'Color','#EDB120','DisplayName','Stretched');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFc),'-', 'Color','#FF00FF','DisplayName','Cut-off');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFg),'-', 'Color','#77AC30','DisplayName','Gamma');
end
ylim(0.5*[-1 1]);
xlabel('TDT Normalized Time (-)'); ylabel('CDF Residual (-)');

% Get just the Basel case for CDF stuff.
ID='Basel';
i=find(strcmpi({Sf.ID},ID));
Sf=Sf(i);

% Plot the Basel CDF.
subplot(122);
% Loop over every case study.
for i=1:length(Sf)

    % Get the value of TDT.
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Omori',Sf(i).Po); Fto=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Exponential',Sf(i).Pe); Fte=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Stretched',Sf(i).Ps); Fts=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Cut-off',Sf(i).Pc); Ftc=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Gamma',Sf(i).Pg); Ftg=Nt/nt;
    Ft=median([Fto Fte Fts Ftc Ftg]);
    
    % Plot the residuals of each model.
    plot(Sf(i).Ts,1:length(Sf(i).Ts),'-k','DisplayName','Data'); hold on;
    plot(Sf(i).t,Sf(i).CDFo,'-', 'Color','#0000FF','DisplayName','Omori');
    plot(Sf(i).t,Sf(i).CDFe,'-', 'Color','#FF0000','DisplayName','Exp');
    plot(Sf(i).t,Sf(i).CDFs,'-', 'Color','#EDB120','DisplayName','Stretched');
    plot(Sf(i).t,Sf(i).CDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off');
    plot(Sf(i).t,Sf(i).CDFg,'-', 'Color','#77AC30','DisplayName','Gamma');
end
xlim([-25 max(xlim())]);
plot(0.2648*31*[1 1],ylim(),'--k'); % Start of flow back [digitized from Häring 2008].
%plot(0.8553*31*[1 1],ylim(),'--k'); % first big decrease in flow back rate.
xlabel('Time (days)'); ylabel('Cumulative Earthquake Counts (-)');



% Subroutines.
function res = residuals1(xd,xf,yf)
  % Local subroutine that calculates the CDF residuals (on fitted x-axis).
  
  yd=interp1(xd+1e-10*(1:length(xd))',1:length(xd),xf,'nearest','extrap');
  res=(yf-yd)/length(xd);
  
end

function res = residuals2(xd,xf,yf)
  % Local subroutine that calculates the CDF residuals (on data x-axis).
  
  %res=xd+1e-10*(1:length(xd))';
  yf=interp1(xf,yf,xd+1e-10*(1:length(xd))','nearest','extrap');
  res=(yf-(1:length(xd))')/length(xd);
  
end