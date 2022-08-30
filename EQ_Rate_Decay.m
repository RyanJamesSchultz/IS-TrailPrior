function [n,N,Nt]=EQ_Rate_Decay(t,type_flag,params)
  % Function that will compute the expected number of earthquake 
  % aftershocks as a function of time, given a model and parameters.
  % See also Schultz et al., [2022].
  % 
  % Reference: 
  % Schultz, R., Ellsworth, W. L., & Beroza, G. C. (2022). Statistical bounds on how induced seismicity stops. Scientific reports, 12(1), 1-11, doi: 10.1038/s41598-022-05216-9.
  %
  % Written by Ryan Schultz.
  
  % Rate Decay Laws.
  no = @(t,K,c,p)     K*(t+c).^(-p);
  ne = @(t,n0,T)      n0*exp(-t/T);
  ns = @(t,Nn,t0,d,q) q*Nn * exp((d/t0).^q) * (1./(t+d)) .* ((t+d)/t0).^q .* exp(-((t+d)/t0).^q);
  nc = @(t,K,c,T,p)   K*(t+c).^(-p).*exp(-t/T);
  ng = @(t,A,lb,la,q) A*t.^-q*gamma(q).*(gammainc(lb*t,q,'lower')-gammainc(la*t,q,'lower'));
  
  % Cumulative Aftershock Counts.
  No = @(t,K,c,p)     (K/(p-1))*(c^(1-p)-(t+c).^(1-p));
  Ne = @(t,n0,T)      T*n0*(1-exp(-t/T));
  Ns = @(t,Nn,t0,d,q) Nn*(1-exp((d/t0).^q)*exp(-((t+d)/t0).^q));
  Nc = @(t,K,c,T,p)   K*T^(1-p)*exp(c/T)*(gammainc(c/T,1-p,'upper')-gammainc((t+c)/T,1-p,'upper'))*gamma(1-p);
  Ng = @(t,A,lb,la,q) (A/(q-1))*((lb^(q-1)-la^(q-1))+la^(q-1)*exp(-t*la)-t.^(1-q).*gammainc(la*t,q,'upper')*gamma(q)-lb^(q-1)*exp(-t*lb)+t.^(1-q).*gammainc(lb*t,q,'upper')*gamma(q));
  
  % Compute the aftershock counts for each of the types.
  if(strcmpi(type_flag,'Omori')) % Modified Omori Power Law [Utsu, 1961].
      K=params(1);
      c=params(2);
      p=params(3);
      n=no(t,K,c,p);
      N=No(t,K,c,p);
      Nt=K*c^(1-p)/(p-1);
  elseif(strcmpi(type_flag,'Exponential')) % Exponential [Burridge & Knopoff, 1967].
      n0=params(1);
      T=params(2);
      n=ne(t,n0,T);
      N=Ne(t,n0,T);
      Nt=T*n0;
  elseif(strcmpi(type_flag,'Stretched')) % Stretched Exponential [Gross & Kisslinger, 1994].
      Nn=params(1);
      t0=params(2);
      d=params(3);
      q=params(4);
      n=ns(t,Nn,t0,d,q);
      N=Ns(t,Nn,t0,d,q);
      Nt=Nn;
      
      % Special handling of the zero-time limit.
      n(t==0)=(q*Nn/d)*(d/t0)^q;
      N(t==0)=0;
  elseif(strcmpi(type_flag,'Cut-off')) % Cut-Off Power Law [Otsuka, 1985].
      K=params(1);
      c=params(2);
      T=params(3);
      p=params(4);
      n=nc(t,K,c,T,p);
      N=Nc(t,K,c,T,p);
      Nt=K*T^(1-p)*exp(c/T)*gammainc(c/T,1-p,'upper')*gamma(1-p);
  elseif(strcmpi(type_flag,'Gamma')) % Gamma Distribution [Narteau et al., 2002].
      A=params(1);
      lb=params(2);
      la=params(3);
      q=params(4);
      n=ng(t,A,lb,la,q);
      N=Ng(t,A,lb,la,q);
      Nt=(A/(q-1))*(lb^(q-1)-la^(q-1));
      
      % Special handling of the zero-time limit.
      n(t==0)=(A/q)*(lb^q-la^q);
      N(t==0)=0;
  end
  
  % Special handling of the infinite-time limit.
  n(t==Inf)=0;
  N(t==Inf)=Nt;
  
end