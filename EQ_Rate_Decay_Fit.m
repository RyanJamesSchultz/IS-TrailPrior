function [pf]=EQ_Rate_Decay_Fit(Td,type_flag,bounds,guess)
  % Function to MLE fit aftershock earthquake times to a given decay model.
  % See also Schultz et al., [2022].
  % 
  % Reference: 
  % Schultz, R., Ellsworth, W. L., & Beroza, G. C. (2022). Statistical bounds on how induced seismicity stops. Scientific reports, 12(1), 1-11, doi: 10.1038/s41598-022-05216-9.
  %
  % Written by Ryan Schultz.
  
  % Ensure that the input times are sorted chronologically.
  Td=Td(:);
  Td=sort(Td);
  
  % Handle the truncation bounds.
  if(isempty(bounds))
      bounds=[0.9*min(Td),1.1*max(Td)];
  end
  
  % Predefine some starting values, define the PDFs & CDFs [Schultz et al., 2022], and then fit the data.
  if(strcmpi(type_flag,'Omori')) % Modified Omori Power Law [Utsu, 1961].
      % Parameters: K, c, p.
      if(isempty(guess))
          p0=[0.01,1.2];
          guess=1;
      else
          p0=guess(2:end);
      end
      % Define the normalized PDF & CDF.
      n = @(t,c,p)     ((t+c).^(-p))/((c.^(1-p))/(p-1));
      N = @(t,c,p)     (1/(p-1))*(c^(1-p)-(t+c).^(1-p))/((c.^(1-p))/(p-1));
      % Fit the normalized PDF/CDF parameters.
      pf=mle(Td,'pdf',n,'start',p0,'cdf',N, 'LowerBound',[1e-10 1+1e-10],'UpperBound',[1e2 1e1],'TruncationBounds',bounds);
      pf=[guess(1),pf];
      % Define the scaled PDF & CDF.
      n = @(t,K,c,p)   K*(t+c).^(-p);
      N = @(t,K,c,p)   (K/(p-1))*(c^(1-p)-(t+c).^(1-p));
      % Define the log-likelihood function [Ogata, 1983; Eqn 6].
      LL = @(t,Tb,K,c,p) sum(log(n(t,K,c,p)))-(N(max(Tb),K,c,p)-N(min(Tb),K,c,p));
      % Fit the scaling parameter.
      pf(1)=fminbnd(@(x) -LL(Td,bounds,x,pf(2),pf(3)), 0,1e20);
      %pf(1)=fminsearch(@(x) -LL(Td,bounds,x,pf(2),pf(3)), guess(1));
      
  elseif(strcmpi(type_flag,'Exponential')) % Exponential [Burridge & Knopoff, 1967].
      % Parameters: n0, T.
      if(isempty(guess))
          p0=[2];
          guess=1;
      else
          p0=guess(2:end);
      end
      % Define the normalized PDF & CDF.
      n = @(t,T)      exp(-t/T)/(T);
      N = @(t,T)      T*(1-exp(-t/T))/(T);
      % Fit the normalized PDF/CDF parameters.
      pf=mle(Td,'pdf',n,'start',p0,'cdf',N, 'LowerBound',[1e-10],'UpperBound',[1e3],'TruncationBounds',bounds);
      pf=[guess(1),pf];
      % Define the scaled PDF & CDF.
      n = @(t,n0,T)   n0*exp(-t/T);
      N = @(t,n0,T)   T*n0*(1-exp(-t/T));
      % Define the log-likelihood function [Ogata, 1983; Eqn 6].
      LL = @(t,Tb,n0,T) sum(log(n(t,n0,T)))-(N(max(Tb),n0,T)-N(min(Tb),n0,T));
      % Fit the scaling parameter.
      pf(1)=fminbnd(@(x) -LL(Td,bounds,x,pf(2)), 0,1e20);
      
  elseif(strcmpi(type_flag,'Stretched')) % Stretched Exponential [Gross & Kisslinger, 1994].
      % Parameters: Nn, t0, d, q.
      if(isempty(guess))
          p0=[150,0.01,0.2];
          guess=1;
      else
          p0=guess(2:end);
      end
      % Define the normalized PDF & CDF.
      n = @(t,t0,d,q) q * exp((d/t0).^q) * (1./(t+d)) .* ((t+d)/t0).^q .* exp(-((t+d)/t0).^q);
      N = @(t,t0,d,q) (1-exp((d/t0).^q)*exp(-((t+d)/t0).^q));
      % Fit the normalized PDF/CDF parameters.
      pf=mle(Td,'pdf',n,'start',p0,'cdf',N, 'LowerBound',[1e-6 0.0 1e-2],'UpperBound',[1e4 1e2 1e1],'TruncationBounds',bounds);
      pf=[guess(1),pf];
      % Define the scaled PDF & CDF.
      n = @(t,Nn,t0,d,q) q*Nn * exp((d/t0).^q) * (1./(t+d)) .* ((t+d)/t0).^q .* exp(-((t+d)/t0).^q);
      N = @(t,Nn,t0,d,q) Nn*(1-exp((d/t0).^q)*exp(-((t+d)/t0).^q));
      % Define the log-likelihood function [Ogata, 1983; Eqn 6].
      LL = @(t,Tb,Nn,t0,d,q) sum(log(n(t,Nn,t0,d,q)))-(N(max(Tb),Nn,t0,d,q)-N(min(Tb),Nn,t0,d,q));
      % Fit the scaling parameter.
      pf(1)=fminbnd(@(x) -LL(Td,bounds,x,pf(2),pf(3),pf(4)), 0,1e20);
      
  elseif(strcmpi(type_flag,'Cut-off')) % Cut-Off Power Law [Otsuka, 1985].
      % Parameters: K, c, T, p.
      if(isempty(guess))
          p0=[0.01,500,0.9];
          guess=1;
      else
          p0=guess(2:end);
      end
      % Define normalized the PDF & CDF.
      n = @(t,c,T,p)     (t+c).^(-p).*exp(-t/T)/(T^(1-p)*exp(c/T)*gammainc(c/T,1-p,'upper')*gamma(1-p));
      N = @(t,c,T,p)     T^(1-p)*exp(c/T)*(gammainc(c/T,1-p,'upper')-gammainc((t+c)/T,1-p,'upper'))*gamma(1-p)/(T^(1-p)*exp(c/T)*gammainc(c/T,1-p,'upper')*gamma(1-p));
      % Fit the normalized PDF/CDF parameters.
      pf=mle(Td,'pdf',n,'start',p0,'cdf',N, 'LowerBound',[1e-15 1e-3 1e-15],'UpperBound',[1e2 1e6 1-1e-10],'TruncationBounds',bounds);
      pf=[guess(1),pf];
      % Define the scaled PDF & CDF.
      n = @(t,K,c,T,p)   K*(t+c).^(-p).*exp(-t/T);
      N = @(t,K,c,T,p)   K*T^(1-p)*exp(c/T)*(gammainc(c/T,1-p,'upper')-gammainc((t+c)/T,1-p,'upper'))*gamma(1-p);
      % Define the log-likelihood function [Ogata, 1983; Eqn 6].
      LL = @(t,Tb,K,c,T,p) sum(log(n(t,K,c,T,p)))-(N(max(Tb),K,c,T,p)-N(min(Tb),K,c,T,p));
      % Fit the scaling parameter.
      pf(1)=fminbnd(@(x) -LL(Td,bounds,x,pf(2),pf(3),pf(4)), 0,1e20);
      
  elseif(strcmpi(type_flag,'Gamma')) % Gamma Distribution [Narteau et al., 2002].
      % Parameters: A, lb, la, q.
      if(isempty(guess))
          p0=[1,1e-6,1.1];
          guess=1;
      else
          p0=guess(2:end);
      end
      % Define the normalized PDF & CDF.
      n = @(t,lb,la,q) t.^-q*gamma(q).*(gammainc(lb*t,q,'lower')-gammainc(la*t,q,'lower'))/((1/(q-1))*(lb^(q-1)-la^(q-1)));
      N = @(t,lb,la,q) (1/(q-1))*((lb^(q-1)-la^(q-1))+la^(q-1)*exp(-t*la)-t.^(1-q).*gammainc(la*t,q,'upper')*gamma(q)-lb^(q-1)*exp(-t*lb)+t.^(1-q).*gammainc(lb*t,q,'upper')*gamma(q))/((1/(q-1))*(lb^(q-1)-la^(q-1)));
      % Fit the normalized PDF/CDF parameters.
      pf=mle(Td,'pdf',n,'start',p0,'cdf',N, 'LowerBound',[p0(1)/2 1e-15 1+1e-10],'UpperBound',[1e3 0.95*p0(1)+1e-10 3],'TruncationBounds',bounds);
      pf=[guess(1),pf];
      pf(3)=min([pf(3) 0.95*pf(2)]);
      % Define the scaled PDF & CDF.
      n = @(t,A,lb,la,q) A*t.^-q*gamma(q).*(gammainc(lb*t,q,'lower')-gammainc(la*t,q,'lower'));
      N = @(t,A,lb,la,q) (A/(q-1))*((lb^(q-1)-la^(q-1))+la^(q-1)*exp(-t*la)-t.^(1-q).*gammainc(la*t,q,'upper')*gamma(q)-lb^(q-1)*exp(-t*lb)+t.^(1-q).*gammainc(lb*t,q,'upper')*gamma(q));
      % Define the log-likelihood function [Ogata, 1983; Eqn 6].
      LL = @(t,Tb,A,lb,la,q) sum(log(n(t,A,lb,la,q)))-(N(max(Tb),A,lb,la,q)-N(min(Tb),A,lb,la,q));
      % Fit the scaling parameter.
      pf(1)=fminbnd(@(x) -LL(Td,bounds,x,pf(2),pf(3),pf(4)), 0,1e20);
      
  end
  
  
end