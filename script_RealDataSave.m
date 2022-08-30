% Script to perform the model fits to all of the real-case data.
clear;

% Load in a predefined structure for convenience.
tableTData;

% Predefine some values.
Nf=1e6;
Ns=5e3;

% Loop over every case considered.
for i=1:length(T)
    
    % Read iteration info to the screen.
    i
    T(i).ID
    
    % Boundary polygons and stuff.
    latB=T(i).latB;
    lonB=T(i).lonB;
    depB=[];
    tB=T(i).Tf;
    mB=[T(i).Mc Inf];
    
    % Get the test sample.
    load(T(i).file);
    
    % Preprocess the data.
    [Ts,cat,bounds]=filtCat(S,latB,lonB,depB,tB,mB);
    t=linspace(0,bounds(2),Ns); t(1)=[];
    
    % Fit each model to the test sample.
    [sf]=Fit_Trailing(Ts,bounds, Nf,t, T(i).Go,T(i).Ge,T(i).Gs,T(i).Gc,T(i).Gg);
    for j=1:5
        [sf]=Fit_Trailing(Ts,bounds, Nf,t, sf.Po,sf.Pe,sf.Ps,sf.Pc,sf.Pg);
    end
    sf.Ts=Ts;
    sf.t=t;
    sf.ID=T(i).ID;

    % Stuff the data into the structure.
    if(i==1)
        Sf=sf;
    else
        Sf(end+1)=sf;
    end

end

% Save the data file.
save('DataFits.mat','Sf','-v7.3');

