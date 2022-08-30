% Script to view the residuals for analysis and QC'ing.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Truncate data?
%i=[17 19 20 21 24 32 37 40 44 46 52 54 55]; % All the 'bad' ones, by KSp (5% confidence).
i=[17 19 20 21 46]; % Halted trail (5%).
i=[8 14 16 22 26]; % Additional 'bad' ones, by KSp (10% confidence).
Sf=Sf(i);

% Plot the CDF residuals.
%figure(1); clf;
%semilogx([1e-3 1e3],[0 0],'-k'); hold on;
% Loop over every case study.
for i=1:length(Sf)
    i
    Sf(i).ID
    W_av=mean([Sf(i).Waic;Sf(i).Wbic]);
    [~,j]=max(W_av);
    if(j==1)
        'Omori'
        yd=Sf(i).CDFo;
    elseif(j==2)
        'Exp'
        yd=Sf(i).CDFe;
    elseif(j==3)
        'Stretch'
        yd=Sf(i).CDFs;
    elseif(j==4)
        'Cut'
        yd=Sf(i).CDFc;
    elseif(j==5)
        'Gamma'
        yd=Sf(i).CDFg;
    end
    length(Sf(i).Ts)
    Sf(i).Waic
    Sf(i).Wbic
    Sf(i).R2b
    Sf(i).KSp
    figure(1); clf;
    semilogx([1e-3 1e3],[0 0],'-k'); hold on;
    %plot(Sf(i).t/max(Sf(i).t),residuals(Sf(i).Ts,Sf(i).t,yd),'-c'); hold on;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Omori',Sf(i).Po); Fto=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Exponential',Sf(i).Pe); Fte=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Stretched',Sf(i).Ps); Fts=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Cut-off',Sf(i).Pc); Ftc=Nt/nt;
    [nt,~,Nt]=EQ_Rate_Decay(1e-8,'Gamma',Sf(i).Pg); Ftg=Nt/nt;
    Ft=median([Fto Fte Fts Ftc Ftg]);
    %Ft=1;
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFo),'-','Color','#0000FF','DisplayName','Omori');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFe),'-','Color','#FF0000','DisplayName','Exp');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFs),'-','Color','#EDB120','DisplayName','Stretched');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFc),'-','Color','#FF00FF','DisplayName','Cut-off');
    semilogx(Sf(i).t/Ft,residuals1(Sf(i).Ts,Sf(i).t,Sf(i).CDFg),'-','Color','#77AC30','DisplayName','Gamma');
    %xlim([0 1]); 
    %ylim(0.5*[-1 1]);
    %xlabel('Sample Number'); ylabel('CDF Residual (-)');
    %legend('Location','northeast');
    pause;
end
%xlim([0 1]);
ylim(0.5*[-1 1]);
xlabel('Normalized Time (-)'); ylabel('CDF Residual (-)');






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