function [R2b,R2r, LL,AICc,BIC,KSp]=Fit_Metrics(te,tf,type_flag,params,bounds)
  % Simple function to compute the perfromance metrics for aftershock
  % distribution data fits.
  %
  % References:
  % Ogata, Y. (1983). Estimation of the parameters in the modified Omori formula for aftershock frequencies by the maximum likelihood procedure. Journal of Physics of the Earth, 31(2), 115-124.
  % Wagenmakers & Farrell (2004). AIC model selection using Akaike weights. Psychonomic bulletin & review, 11(1), 192-196, doi: 10.3758/BF03206482.
  %
  % Written by Ryan Schultz.
  
  % Sort the input data.
  te=sort(te);
  tf=sort(tf);
  
  % Define the number of parameters for each model.
  n=length(te);
  if(strcmpi(type_flag,'Omori')) % Modified Omori Power Law [Utsu, 1961].
      K=1+3;
  elseif(strcmpi(type_flag,'Exponential')) % Exponential [Burridge & Knopoff, 1967].
      K=1+2;
  elseif(strcmpi(type_flag,'Stretched')) % Stretched Exponential [Gross & Kisslinger, 1994].
      K=1+4;
  elseif(strcmpi(type_flag,'Cut-off')) % Cut-Off Power Law [Otsuka, 1985].
      K=1+4;
  elseif(strcmpi(type_flag,'Gamma')) % Gamma Distribution [Narteau et al., 2002].
      K=1+4;
  end
  
  % Compute the empirical PDFs & CDFs.
  [PDFe,CDFe,tb]=Get_Empirical_DF(te,0);
  [PDFf,CDFf,~]=Get_Empirical_DF(tf,tb);
  
  % Compute the adjusted goodness-of-fit (R²) statistic (binned).
  w=sqrt(PDFe); w=w/sum(w);
  Ye=PDFe;
  Yf=PDFf;
  Yb=sum(w.*Ye);
  SStot=sum(w.*((Ye-Yb).^2));
  SSres=sum(w.*((Ye-Yf).^2));
  R2b=1-(SSres/SStot);
  R2b=1-((1-R2b)*(n-1)/(n-K-1));
  
  % Compute the adjusted goodness-of-fit (R²) statistic (unbinned).
  Ye=1./diff(te);
  [Yf,~,~]=EQ_Rate_Decay(te,type_flag,params); Yf=Yf(2:end);
  I=~isinf(Ye); Ye=Ye(I); Yf=Yf(I);
  w=ones(size(Ye)); w=w/sum(w);
  Yb=sum(w.*Ye);
  SStot=sum(w.*((Ye-Yb).^2));
  SSres=sum(w.*((Ye-Yf).^2));
  R2r=1-(SSres/SStot);
  R2r=1-((1-R2r)*(n-1)/(n-K-1));
  
  % Compute the log-likelihood [Ogata, 1983].
  [PDFp,~,Np]=EQ_Rate_Decay(te,type_flag,params);
  [~,CDFp,~]=EQ_Rate_Decay(bounds,type_flag,params);
  LL=sum(log(PDFp))-(CDFp(2)-CDFp(1)); % Eqn 6.
  
  % Compute the AIC & BIC statistics [Wagenmakers & Farrell, 2004].
  %AIC=2*K-2*LL;
  AICc=2*K+(2*K*(K+1)/(n-K-1))-2*LL;
  BIC=K*log(n)-2*LL;
  
  % Kolmogorov-Smirnov test p-value.
  [~,KSp]=kstest2(te,tf);
  
  % In case there's too few samples for the number of fit parameters.
  if(K>=n)
      R2b=-Inf;
      R2r=-Inf;
      AICc=Inf;
      BIC=Inf;
  end

end



% Subroutine
function [PDF,CDF,tb]=Get_Empirical_DF(te,tb)
  % Simple subrountine to compute the empirical PDF and CDFs of data.
  %
  % Written by Ryan Schultz.
  
  % Compute the empirical PDF & CDF.
  te=sort(te);
  CDF=(1:length(te))/length(te);
  if(tb==0)
      [PDF,tb]=histcounts(te,round(2*sqrt(length(te))),'Normalization','pdf');
  else
      [PDF,tb]=histcounts(te,tb,'Normalization','pdf');
  end
  
  
end