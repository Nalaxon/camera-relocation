%Autor: Thomas Pitsch, Bernhard Rapp
%Date:  23.10.2014
%Version 0.1
%Use Piotr's Matlab Toolbox [1] to train a regression random forest for camera
%relocation
%[1] http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

%do not forget load toolbox ;-)

%define forst configuration
N=1000; d=2; H=5;
[xs0, hs0, xs1, hs1] = demoGenData(N, N,H,d,1,1);
xs0=single(xs0); xs1=single(xs1);

%train forest
pTrain={'maxDepth', 50, 'F1', 2 'M', 150, 'minChild', 5, 'H', H};

%train forst
tic, forest=forestTrain(xs0, hs0, pTrain{:}); toc

%apply forst on first dataset
hsPr0 = forestApply(xs0, forest);

%apply forst on second dataset
hsPr1 = forestApply(xs1, forest);

%some output they thought is usefull
e0=mean(hsPr0~=hs0); e1=mean(hsPr1~=hs1);
fprintf('errors trn=%f tst=%f\n', e0, e1); figure(1);
figure(1);
subplot(2,2,1); visualizeData(xs0,2,hs0);
subplot(2,2,2); visualizeData(xs0,2,hsPr0);
subplot(2,2,3); visualizeData(xs1,2,hs1);
subplot(2,2,4); visualizeData(xs1,2,hsPr1);

