function [Sc]=Fit_Trailing2(Ti,bounds,Tc,Nm,Tf,guess_o,guess_e,guess_s,guess_c,guess_g)
  % Simple wrapper to fit all of the aftershock models to a trailing
  % seismicity dataset in a two-part piecwise fashion.
  %
  % Please forgive the crudeness of this routine, I hacked it together just
  % to quickly answer a question.
  % 
  % Written by Ryan Schultz.
  
  % Cut the data into two parts.
  Ti1=Ti(Ti<=Tc);
  Ti2=Ti(Ti>Tc);
  
  % Fit each of the two peices.
  [S1]=Fit_Trailing(Ti1,[bounds(1) Tc], Nm,Tf, guess_o,guess_e,guess_s,guess_c,guess_g);
  [S2]=Fit_Trailing(Ti2,[Tc bounds(2)], Nm,Tf, guess_o,guess_e,guess_s,guess_c,guess_g);
  Sc.S1=S1;
  Sc.S2=S2;
  
  % Make the composite PDFs/CDFs.
  Sc.PDFo=[S1.PDFo(Tf<=Tc),S2.PDFo(Tf>Tc)];
  Sc.PDFe=[S1.PDFe(Tf<=Tc),S2.PDFe(Tf>Tc)];
  Sc.PDFs=[S1.PDFs(Tf<=Tc),S2.PDFs(Tf>Tc)];
  Sc.PDFc=[S1.PDFc(Tf<=Tc),S2.PDFc(Tf>Tc)];
  Sc.PDFg=[S1.PDFg(Tf<=Tc),S2.PDFg(Tf>Tc)];
  
  % Make the composite PDFs/CDFs.
  j=find(Tf>Tc,1,'first');
  Sc.CDFo=[S1.CDFo(Tf<=Tc),S2.CDFo(Tf>Tc)];  Sc.CDFo(j:end)=Sc.CDFo(j:end)+(S1.CDFo(j)-Sc.CDFo(j));
  Sc.CDFe=[S1.CDFe(Tf<=Tc),S2.CDFe(Tf>Tc)];  Sc.CDFe(j:end)=Sc.CDFe(j:end)+(S1.CDFe(j)-Sc.CDFe(j));
  Sc.CDFs=[S1.CDFs(Tf<=Tc),S2.CDFs(Tf>Tc)];  Sc.CDFs(j:end)=Sc.CDFs(j:end)+(S1.CDFs(j)-Sc.CDFs(j));
  Sc.CDFc=[S1.CDFc(Tf<=Tc),S2.CDFc(Tf>Tc)];  Sc.CDFc(j:end)=Sc.CDFc(j:end)+(S1.CDFc(j)-Sc.CDFc(j));
  Sc.CDFg=[S1.CDFg(Tf<=Tc),S2.CDFg(Tf>Tc)];  Sc.CDFg(j:end)=Sc.CDFg(j:end)+(S1.CDFg(j)-Sc.CDFg(j));
  
  % Get the composite log-likelihoods.
  Sc.LL(1)=sum(log(interp1(Tf,Sc.PDFo,Ti,'linear')))-(interp1([0,Tf],[0,Sc.CDFo],bounds(1))-interp1(Tf,Sc.CDFo,bounds(2)));
  Sc.LL(2)=sum(log(interp1(Tf,Sc.PDFe,Ti,'linear')))-(interp1([0,Tf],[0,Sc.CDFe],bounds(1))-interp1(Tf,Sc.CDFe,bounds(2)));
  Sc.LL(3)=sum(log(interp1(Tf,Sc.PDFs,Ti,'linear')))-(interp1([0,Tf],[0,Sc.CDFs],bounds(1))-interp1(Tf,Sc.CDFs,bounds(2)));
  Sc.LL(4)=sum(log(interp1(Tf,Sc.PDFc,Ti,'linear')))-(interp1([0,Tf],[0,Sc.CDFc],bounds(1))-interp1(Tf,Sc.CDFc,bounds(2)));
  Sc.LL(5)=sum(log(interp1(Tf,Sc.PDFg,Ti,'linear')))-(interp1([0,Tf],[0,Sc.CDFg],bounds(1))-interp1(Tf,Sc.CDFg,bounds(2)));
  
  % Get the AIC & BIC scores.
  K=[2*3+1+1, 2*2+1+1, 2*4+1+1, 2*4+1+1, 2*4+1+1];
  n=length(Ti);
  Sc.AIC=2*K+(2*K.*(K+1)./(n-K-1))-2*Sc.LL;
  Sc.BIC=K*log(n)-2*Sc.LL;

  % Compute the delta AIC and BIC values.
  Sc.dAIC=Sc.AIC-min(Sc.AIC);
  Sc.dBIC=Sc.BIC-min(Sc.BIC);
  
  % Compute model weights, based on AIC and BIC.
  Sc.Waic=exp(-Sc.dAIC/2)/sum(exp(-Sc.dAIC/2));
  Sc.Wbic=exp(-Sc.dBIC/2)/sum(exp(-Sc.dBIC/2));
  
  % Get the Kolmogorov-Smirnov test p-value.
  Tm=EQ_Rate_Decay_Rand2([0,Tf],[0,Sc.CDFo],[Nm 1]); [~,KSp]=kstest2(Ti,Tm); Sc.KSp(1)=log10(KSp);
  Tm=EQ_Rate_Decay_Rand2([0,Tf],[0,Sc.CDFe],[Nm 1]); [~,KSp]=kstest2(Ti,Tm); Sc.KSp(2)=log10(KSp);
  Tm=EQ_Rate_Decay_Rand2([0,Tf],[0,Sc.CDFs],[Nm 1]); [~,KSp]=kstest2(Ti,Tm); Sc.KSp(3)=log10(KSp);
  Tm=EQ_Rate_Decay_Rand2([0,Tf],[0,Sc.CDFc],[Nm 1]); [~,KSp]=kstest2(Ti,Tm); Sc.KSp(4)=log10(KSp);
  Tm=EQ_Rate_Decay_Rand2([0,Tf],[0,Sc.CDFg],[Nm 1]); [~,KSp]=kstest2(Ti,Tm); Sc.KSp(5)=log10(KSp);
  
end




function [T]=EQ_Rate_Decay_Rand2(t,N,Nr)
  % Subroutine that randomly draws earthquake times from a given rate decay model.
  % 
  % Written by Ryan Schultz.
  
  % Get the time horizon and the cumulative number of events.
  %[~,N,~]=EQ_Rate_Decay(t,type_flag,params);
  
  % Normalize to a CDF.
  N=N-N(1);
  N=N/N(end);
  
  % Make sure they're sorted and unique.
  [~,I,~]=unique(N);
  N=N(I);
  t=t(I);
  
  % Get the random times.
  r=rand(Nr);
  T=interp1(N,t,r,'linear');
  
  % Sort the output aftershock times.
  T=sort(T);
  
end