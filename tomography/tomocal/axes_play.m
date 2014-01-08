% Try to unscrambles axes
%% generate some random axes.
ar=0.1; % range of angle errors allowed
sr=0.1; % range of losses allowed.
t1 = rand(1)*ar - ar/2;
t2 = rand(1)*ar - ar/2;
t1=0; t2=0; % Assume we can calibrate this error out from singlet-data.
p1 = (rand(1)*ar - ar/2);
p2 = (rand(1)*ar - ar/2)*0;
se1 = 1-sr*rand(1);
se2 = 1-sr*rand(1);
axes=[ 1 0 0 ; ...
       sin(t1) cos(t1)*cos(p1) cos(t1)*sin(p1) ; ...
       sin(t1) cos(t2)*sin(p2) cos(t2)*cos(p2)];
axes(2,:)=axes(2,:) * se1;
axes(3,:)=axes(3,:) * se2;
axes=axes';

%% Generate some data.
% jprime, including weird state prep errors..
jprime=generate_rotation(rmat([0 1 0],0.05)*[0;1;0], [1;0;1], 8*pi, 12, 80);
% a mostly-exchange based rotation
j=generate_rotation(rmat([1 0 0],0.05)*[0;1;0], [1;0;0.05], 8*pi, 10, 60);
%dbz
dbz=generate_rotation(rmat([1 0 0],0.03)*[1;0;0], [0.05;0;1], 8*pi, 8, 32*4);
jprime=skew_axes(jprime,axes);
j=skew_axes(j,axes);
dbz=skew_axes(dbz,axes);

%% Add some noise
j.data=j.data+randn(size(j.data))*1e-2;
jprime.data=jprime.data+randn(size(jprime.data))*1e-2;
dbz.data=dbz.data+randn(size(dbz.data))*1e-2;

%% Add some sensor biases
bias=[0.01 -0.01 0.01]*1;
j.data=j.data+repmat(bias,size(j.data,1),1);
jprime.data=jprime.data+repmat(bias,size(jprime.data,1),1);
%%
% Try to fit the data
mfdata(1).x=jprime.t;
mfdata(1).y=jprime.data(:);  % We'll need to repack this is the fit function
[mfmodel(1).yfn mfmodel(1).fn] = skew_fit_function([5 6]);

mfdata(2).x=j.t;
mfdata(2).y=j.data(:);  % We'll need to repack this is the fit function
[mfmodel(2).yfn mfmodel(2).fn] = skew_fit_function([7 8]);

mfdata(3).x=dbz.t;
mfdata(3).y=dbz.data(:);  % We'll need to repack this is the fit function
[mfmodel(3).yfn mfmodel(3).fn] = skew_fit_function([9 10]);

p0=[ .95 1e-3 0 .95  1 10 1 10 1 10];
pmask=[1 1 0 1 1 1 1 1 1 1];

[p, chisq] = mfitwrap(mfdata,mfmodel,p0,[],[0 0 0 0 1 1 1 1 1 1] & pmask);
pause
%[p, chisq] = mfitwrap(mfdata,mfmodel,p,[],[1 1 0 1 0 0 0 0 0 0] & pmask);
%pause
[p, chisq] = mfitwrap(mfdata,mfmodel,p ,[], pmask);
pm=zeros(3);
pm(2:3,2:3)=reshape(p(1:4),2,2);
pm(1,1)=1;
pm*axes
pm
