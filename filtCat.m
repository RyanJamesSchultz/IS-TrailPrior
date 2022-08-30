function [Ts,cat,bounds,dT]=filtCat(S,lat_L,lon_L,dep_L,T_L,M_L)
  % Simple function to filter and pre-process a catalogue structure.
  
  % Pull out the catalogue.
  cat=S.cat;

  % Filter spatially (lateral).
  if(~isempty(lat_L))
      I=inpolygon(cat.lon,cat.lat,lon_L,lat_L);
      cat.time=cat.time(I);
      cat.mag=cat.mag(I);
      cat.lat=cat.lat(I);
      cat.lon=cat.lon(I);
      cat.dep=cat.dep(I);
  end
  
  % Filter spatially (depth).
  if(~isempty(dep_L))
      I=(cat.dep>=min(dep_L))&(cat.dep<=max(dep_L));
      cat.time=cat.time(I);
      cat.mag=cat.mag(I);
      cat.lat=cat.lat(I);
      cat.lon=cat.lon(I);
      cat.dep=cat.dep(I);
  end
      
  % Filter temporally.
  if(~isempty(T_L))
      I=(cat.time>min(T_L))&(cat.time<max(T_L));
      cat.time=cat.time(I);
      cat.mag=cat.mag(I);
      cat.lat=cat.lat(I);
      cat.lon=cat.lon(I);
      cat.dep=cat.dep(I);
  end
  
  % Filter by magnitudes.
  if(~isempty(M_L))
      I=(cat.mag>=min(M_L))&(cat.mag<=max(M_L));
      cat.time=cat.time(I);
      cat.mag=cat.mag(I);
      cat.lat=cat.lat(I);
      cat.lon=cat.lon(I);
      cat.dep=cat.dep(I);
  end
  
  % Sort chronologically.
  [~,I]=sort(cat.time);
  cat.time=cat.time(I);
  cat.mag=cat.mag(I);
  cat.lat=cat.lat(I);
  cat.lon=cat.lon(I);
  cat.dep=cat.dep(I);
  
  % Get sample times relative from the given start time.
  Ts=cat.time-min(T_L); % Units of days.
  
  % Get the sampling rate of the injection data.
  dT=diff(S.inj.time);
  
  % Define the time sample lower bound.
  bounds(1)=0.9*min(Ts);
  if(isfinite(min(T_L)))
      bounds(1)=min([bounds(1) 1e-5]);
  end
  
  % Define the time sample upper bound.
  bounds(2)=1.1*max(Ts);
  if(isfinite(max(T_L)))
      bounds(2)=min([bounds(2) max(T_L)-min(T_L)]);
  end
  
return


