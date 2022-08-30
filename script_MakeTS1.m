% Script to start making Table S1.
clear;

% Load the in the data structures.
tableTData;
load('DataFits.mat','Sf');

% Define some variables.
n=length(Sf);

% Loop over all of the entries.
for i=1:n
    %Sf(i).ID
    
    % Find the table entry that matches.
    j=find(strcmpi({T.ID},Sf(i).ID));
    
    % Get some start/stop time info.
    Ts=datestr(min(T(j).Tf),'yyyy-mm-ddTHH:MM:SS.FFF');
    if(isfinite(max(T(j).Tf)))
        Te=datestr(max(T(j).Tf),'yyyy-mm-ddTHH:MM:SS.FFF');
    else
        Te='3000-01-01T01:01:01.000';
    end
    
    % Print relevant data:
    
    % Pre-processing stuff.
    % ID, Grade, Mc, Ntrail, S, T,
    fprintf(1,'%s, %s, %0.2f, %d, ',Sf(i).ID, T(j).ID, T(j).Mc, length(Sf(i).Ts) );
    % S, T,
    fprintf(1,'%s, %s, ',Ts,Te );
    
    % Parameter fits.
    % Omori params: K, c, p.
    fprintf(1,'%0.4e, %0.4e, %0.4e, ',Sf(i).Po(1), Sf(i).Po(2), Sf(i).Po(3) );
    % Exponential params: n0, T.
    fprintf(1,'%0.4e, %0.4e, ',Sf(i).Pe(1), Sf(i).Pe(2) );
    % Stretch params: Nn, t0, d, q.
    fprintf(1,'%0.4e, %0.4e, %0.4e, %0.4e, ',Sf(i).Ps(1), Sf(i).Ps(2), Sf(i).Ps(3), Sf(i).Ps(4) );
    % Cut-off params: K, c, T, p.
    fprintf(1,'%0.4e, %0.4e, %0.4e, %0.4e, ',Sf(i).Pc(1), Sf(i).Pc(2), Sf(i).Pc(3), Sf(i).Pc(4) );
    % Gamma params: A, lb, la, q.
    fprintf(1,'%0.4e, %0.4e, %0.4e, %0.4e, ',Sf(i).Pg(1), Sf(i).Pg(2), Sf(i).Pg(3), Sf(i).Pg(4) );
    
    % Performance metrics (Ordering: Omo,Exp,Str,Cut,Gam).
    % AIC, WAIC, BIC, WBIC scores.
    fprintf(1,'%0.4e, %0.4e, %0.4e, %0.4e, %0.4e, ',Sf(i).AIC(1),  Sf(i).AIC(2),  Sf(i).AIC(3),  Sf(i).AIC(4),  Sf(i).AIC(5)  );
    fprintf(1,'%0.4f, %0.4f, %0.4f, %0.4f, %0.4f, ',Sf(i).Waic(1), Sf(i).Waic(2), Sf(i).Waic(3), Sf(i).Waic(4), Sf(i).Waic(5) );
    fprintf(1,'%0.4e, %0.4e, %0.4e, %0.4e, %0.4e, ',Sf(i).BIC(1),  Sf(i).BIC(2),  Sf(i).BIC(3),  Sf(i).BIC(4),  Sf(i).BIC(5)  );
    fprintf(1,'%0.4f, %0.4f, %0.4f, %0.4f, %0.4f, ',Sf(i).Wbic(1), Sf(i).Wbic(2), Sf(i).Wbic(3), Sf(i).Wbic(4), Sf(i).Wbic(5) );
    % log10(KSp) values.
    fprintf(1,'%0.4f, %0.4f, %0.4f, %0.4f, %0.4f, ',Sf(i).KSp(1),  Sf(i).KSp(2),  Sf(i).KSp(3),  Sf(i).KSp(4),  Sf(i).KSp(5) );
    % R2 values.
    fprintf(1,'%0.4f, %0.4f, %0.4f, %0.4f, %0.4f, ',Sf(i).R2b(1),  Sf(i).R2b(2),  Sf(i).R2b(3),  Sf(i).R2b(4),  Sf(i).R2b(5) );
    
    % End line.
    fprintf(1,'\n');
    
end

