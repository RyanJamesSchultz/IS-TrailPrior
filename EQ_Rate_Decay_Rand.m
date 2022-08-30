function [T]=EQ_Rate_Decay_Rand(t,Nr,type_flag,params)
  % Function that randomly draws earthquake times from a given rate decay model.
  % 
  % Written by Ryan Schultz.
  
  % Get the time horizon and the cumulative number of events.
  [~,N,~]=EQ_Rate_Decay(t,type_flag,params);
  
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