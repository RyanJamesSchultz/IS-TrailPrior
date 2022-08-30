function fxnTimes(h,Xlim)
  % Stupid fxn to get end time of a current injection up to the start time 
  % of the next injection in the current view.  Used to get S & T times.
  %
  % Written by Ryan Schultz.
  
  % Load in the data from the current view.
  x=h(1).XData;
  y=h(1).YData;
  
  % Keep on the relevant portions.
  I=(x>=min(Xlim))&(x<=max(Xlim));
  x=x(I);
  y=y(I);
  %x(1)
  %datestr(x(1))

  % Get the start/end times.
  x1=x(find(y<=0, 1, 'first' ));
  x2=x(find(y<=0, 1, 'last'  ));
  
  % Print info to the screen.
  fprintf(1,'%12.5f %12.5f\n',x1,x2);
  datestr([x1 x2])
  
end