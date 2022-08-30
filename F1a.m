% Script to make Figure 1a.
clear;

% Predefine some values.
dt=1e-3;
t=0:dt:1e3; t(1)=[];
Nsamp=1e2;

% Omori.
c1=0.1; p1=1.001;
K=10*Nsamp/(c1^-p1);
params=[K,c1,p1];
[PDFo,CDFo,No]=EQ_Rate_Decay(t,'Omori',params);

% Exponential.
n0=6*Nsamp; tau=50;
params=[n0,tau];
[PDFe,CDFe,Ne]=EQ_Rate_Decay(t,'Exponential',params);

% Stretched exponential.
d=0.06; t0=150; q1=0.01;
Ns=3*Nsamp/((q1/d)*(d/t0)^q1); 
params=[Ns,t0,d,q1];
[PDFs,CDFs,Ns]=EQ_Rate_Decay(t,'Stretched',params);

% Cut-off power law.
c2=0.03; T=100; p2=0.999;
K=1.75*Nsamp/(c2^-p2);
params=[K,c2,T,p2];
[PDFc,CDFc,Nc]=EQ_Rate_Decay(t,'Cut-off',params);

% Gamma.
lb=100; la=1/30; q2=1.001;
A=Nsamp/((1/q2)*(lb^q2-la^q2));
params=[A,lb,la,q2];
[PDFg,CDFg,Ng]=EQ_Rate_Decay(t,'Gamma',params);

figure(1); clf;
%subplot(211);
loglog(t,PDFo,'-', 'Color','#0000FF','DisplayName','Modified-Omori'); hold on;
loglog(c1,interp1(t,PDFo,c1,'linear'),'o','MarkerEdgeColor','#0000FF');
loglog(t,PDFe,'-', 'Color','#FF0000','DisplayName','Exponential');
loglog(tau,interp1(t,PDFe,tau,'linear'),'o','MarkerEdgeColor','#FF0000');
loglog(t,PDFs,'-', 'Color','#EDB120','DisplayName','Stretched Exponential');
loglog([d t0],interp1(t,PDFs,[d t0],'linear'),'o','MarkerEdgeColor','#EDB120');
loglog(t,PDFc,'-', 'Color','#FF00FF','DisplayName','Cut-off Power Law');
loglog([c2 T],interp1(t,PDFc,[c2 T],'linear'),'o','MarkerEdgeColor','#FF00FF');
loglog(t,PDFg,'-', 'Color','#77AC30','DisplayName','Gamma');
loglog([1/lb 1/la],interp1(t,PDFg,[1/lb 1/la],'linear'),'o','MarkerEdgeColor','#77AC30');
xlabel('Time (days)'); ylabel('Trailing Event Rate (events/day)');
legend('Location','southwest');
ylim([1e-3 1.5*max(PDFo)]);
%subplot(212);
%loglog(t,CDFo,'-b','DisplayName','Modified-Omori'); hold on;
%loglog(t,CDFe,'-r','DisplayName','Exponential');
%loglog(t,CDFs,'-y','DisplayName','Stretched Exponential');
%loglog(t,CDFc,'-m','DisplayName','Cut-off');
%loglog(t,CDFg,'-g','DisplayName','Gamma');
%xlabel('Time (days)'); ylabel('Cumulative Counts');
%legend('Location','southeast');
%ylim([0.7 1.5*Nsamp]);
