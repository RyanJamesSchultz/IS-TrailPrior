function [S]=Fit_Trailing(Ti,bounds,Nm,Tf,guess_o,guess_e,guess_s,guess_c,guess_g)
  % Simple wrapper to fit all of the aftershock models to a trailing
  % seismicity dataset and return a structure to house all of the metrics.
  % Performance metric vectors are ordered: Omori, Exp, Stretch, Cut, Gamma.
  % 
  % Written by Ryan Schultz.
  
  % Sort the input data.
  Ti=sort(Ti);
  
  % Predefine the data structure.
  S=struct('R2b',[],'R2r',[],'LL',[],'AIC',[],'BIC',[],'dAIC',[],'dBIC',[],'Waic',[],'Wbic',[],'KSp',[],'Po',[],'PDFo',[],'CDFo',[],'Pe',[],'PDFe',[],'CDFe',[],'Ps',[],'PDFs',[],'CDFs',[],'Pc',[],'PDFc',[],'CDFc',[],'Pg',[],'PDFg',[],'CDFg',[]);
  
  % Predefine initial guesses, if nothing is given.
  if(isempty(guess_o))
      guess_o=[1,0.01,1.2];
  end
  if(isempty(guess_e))
      guess_e=[1,3];
  end
  if(isempty(guess_s))
      guess_s=[1,0.01,0.005,0.15];
  end
  if(isempty(guess_c))
      guess_c=[1,0.001,450,0.99];
  end
  if(isempty(guess_g))
      guess_g=[1,1,0.001,1.1];
  end
  
  %%% The Modified Omori Power Law [Utsu, 1961].
  % Fit the model.
  [params]=EQ_Rate_Decay_Fit(Ti,'Omori',bounds,guess_o);
  [PDFo,CDFo,No]=EQ_Rate_Decay(Tf,'Omori',params);
  % Get the performance metrics.
  Tm=EQ_Rate_Decay_Rand(Tf,[Nm 1],'Omori',params);
  [R2b,R2r,LL,AICc,BIC,KSp]=Fit_Metrics(Ti,Tm,'Omori',params,bounds);
  % Save the outputs in the data structure.
  S.R2b(1)=R2b; S.R2r(1)=R2r; S.LL(1)=LL; S.AIC(1)=AICc; S.BIC(1)=BIC; S.KSp(1)=log10(KSp);
  S.PDFo=PDFo; S.CDFo=CDFo;
  S.Po=params;
  
  %%% The Exponential Decay [Burridge & Knopoff, 1967].
  % Fit the model.
  [params]=EQ_Rate_Decay_Fit(Ti,'Exponential',bounds,guess_e);
  [PDFe,CDFe,Ne]=EQ_Rate_Decay(Tf,'Exponential',params);
  % Get the performance metrics.
  Tm=EQ_Rate_Decay_Rand(Tf,[Nm 1],'Exponential',params);
  [R2b,R2r,LL,AICc,BIC,KSp]=Fit_Metrics(Ti,Tm,'Exponential',params,bounds);
  % Save the outputs in the data structure.
  S.R2b(2)=R2b; S.R2r(2)=R2r; S.LL(2)=LL; S.AIC(2)=AICc; S.BIC(2)=BIC; S.KSp(2)=log10(KSp);
  S.PDFe=PDFe; S.CDFe=CDFe;
  S.Pe=params;
  
  %%% The Stretched Exponential Decay [Gross & Kisslinger, 1994].
  % Fit the model.
  [params]=EQ_Rate_Decay_Fit(Ti,'Stretched',bounds,guess_s);
  [PDFs,CDFs,Ns]=EQ_Rate_Decay(Tf,'Stretched',params);
  % Get the performance metrics.
  Tm=EQ_Rate_Decay_Rand(Tf,[Nm 1],'Stretched',params);
  [R2b,R2r,LL,AICc,BIC,KSp]=Fit_Metrics(Ti,Tm,'Stretched',params,bounds);
  % Save the outputs in the data structure.
  S.R2b(3)=R2b; S.R2r(3)=R2r; S.LL(3)=LL; S.AIC(3)=AICc; S.BIC(3)=BIC; S.KSp(3)=log10(KSp);
  S.PDFs=PDFs; S.CDFs=CDFs;
  S.Ps=params;
  
  %%% The Cut-Off Power Law [Otsuka, 1985].
  % Fit the model.
  [params]=EQ_Rate_Decay_Fit(Ti,'Cut-off',bounds,guess_c);
  [PDFc,CDFc,Nc]=EQ_Rate_Decay(Tf,'Cut-off',params);
  % Get the performance metrics.
  Tm=EQ_Rate_Decay_Rand(Tf,[Nm 1],'Cut-off',params);
  [R2b,R2r,LL,AICc,BIC,KSp]=Fit_Metrics(Ti,Tm,'Cut-off',params,bounds);
  % Save the outputs in the data structure.
  S.R2b(4)=R2b; S.R2r(4)=R2r; S.LL(4)=LL; S.AIC(4)=AICc; S.BIC(4)=BIC; S.KSp(4)=log10(KSp);
  S.PDFc=PDFc; S.CDFc=CDFc;
  S.Pc=params;
  
  %%% The Gamma Distribution [Narteau et al., 2002].
  % Fit the model.
  [params]=EQ_Rate_Decay_Fit(Ti,'Gamma',bounds,guess_g);
  [PDFg,CDFg,Ng]=EQ_Rate_Decay(Tf,'Gamma',params);
  % Get the performance metrics.
  Tm=EQ_Rate_Decay_Rand(Tf,[Nm 1],'Gamma',params);
  [R2b,R2r,LL,AICc,BIC,KSp]=Fit_Metrics(Ti,Tm,'Gamma',params,bounds);
  % Save the outputs in the data structure.
  S.R2b(5)=R2b; S.R2r(5)=R2r; S.LL(5)=LL; S.AIC(5)=AICc; S.BIC(5)=BIC; S.KSp(5)=log10(KSp);
  S.PDFg=PDFg; S.CDFg=CDFg;
  S.Pg=params;
  
  % Compute the delta AIC and BIC values.
  S.dAIC=S.AIC-min(S.AIC);
  S.dBIC=S.BIC-min(S.BIC);
  
  % Compute model weights, based on AIC and BIC.
  S.Waic=exp(-S.dAIC/2)/sum(exp(-S.dAIC/2));
  S.Wbic=exp(-S.dBIC/2)/sum(exp(-S.dBIC/2));
  
end