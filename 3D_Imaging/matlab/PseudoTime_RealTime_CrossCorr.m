cd ~/Projects/PMO/MeasuringFitnessPerClone/code/3D_Imaging/matlab/
addpath cxcorr
FoV='FoFX007_220523_brightfield';
csvf=['~/Projects/PMO/MeasuringFitnessPerClone/data/GastricCancerCL/3Dbrightfield/NCI-N87/A08_Pseudotime/',FoV,'.csv'];

%% Read data
csv=readtable(csvf);
a=table2array(csv(:,'time_since_division'));
b=table2array(csv(:,'pseudotime'));
% a=[1,2,3,4,5,6];
% b=[7,8,3,4,5,6];

[x,c]=cxcorr(a',b');
c=max(c);
%% compare to cross correlation obtained from all random permuations of a
A_=perms(a)';
for a_=A_
    [~,c_]=cxcorr(a_',b');
    c=[c,max(c_)];
end
[~,p]=ttest2(c(1),c(2:length(c)));

[~,i]=max(c);
i=i-2; %% 0 shift possible

%% Shift data to maximize correlation
b1=b((length(b)-i):length(b));
b2=b(1:(length(b)-i-1));

%% Append to end of matrix:
csv.(size(csv, 2)+1) = [b1;b2];
% Change column name if needed
csv.Properties.VariableNames{size(csv, 2)} = 'pseudotime_shifted';
writetable(csv,strrep(csvf,'.csv','_matlabOut.csv'))
