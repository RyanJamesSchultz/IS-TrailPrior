% Script to assist in finding best fit cases.
clear;

% Load the in the data structure.
load('DataFits.mat','Sf');

% Index for model choice (1 omo, 2 exp, 3 str, 4 cut, 5 gam).
i=4;

% Keep only good data.
If=arrayfun(@(S) S.KSp(i),Sf)>=log10(0.05);
Sf=Sf(If);

% Get the ordering of best to worst, for each metric.
[~,Ia]=sort(arrayfun(@(S) S.Waic(i),Sf),'descend');
[~,Ib]=sort(arrayfun(@(S) S.Wbic(i),Sf),'descend');
[~,Ic]=sort(arrayfun(@(S) (S.Waic(i)+S.Wbic(i))/2,Sf),'descend');
[~,Ir]=sort(arrayfun(@(S) S.R2b(i),Sf),'descend');
[~,Ik]=sort(arrayfun(@(S) S.KSp(i),Sf),'descend');

% List rankings.
Ia(1:5)
Ib(1:5)
Ic(1:5)
Ir(1:5)
Ik(1:5)

% Give the top three from the AIC/BIC compromise.
Sf(Ic(1))
Sf(Ic(2))
Sf(Ic(3))