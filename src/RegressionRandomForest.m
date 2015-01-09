%Autor: Thomas Pietsch, Bernhard Rapp
%Date:  23.10.2014
%Version 0.1
%Use Piotr's Matlab Toolbox [1] to train a regression random forest for camera
%relocation
%[1] http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

%do not forget load toolbox ;-)



%toolboxCompile;


%define forst configuration
N=10000; d=5; H=5;
[xs0, hs0, xs1, hs1] = demoGenData(N, N,H,d,1,1);
xs0=single(xs0); xs1=single(xs1);

%train forest
pTrain={'maxDepth', 50, 'F1', 2 'M', 1500, 'minChild', 1, 'H', H, 'split', 'custom'};

%train forst
forest=forestRegTrain(xs0, hs0, pTrain{:});

%apply forst on first dataset
hsPr0 = forestRegApply(xs0, forest);

%apply forst on second dataset
hsPr1 = forestRegApply(xs1, forest);

%some output they thought is usefull
e0=mean(hsPr0~=hs0); e1=mean(hsPr1~=hs1);
fprintf('errors trn=%f tst=%f\n', e0, e1); figure(1);
figure(1);
subplot(3,2,1); visualizeData(xs0,2,hs0);
subplot(3,2,2); visualizeData(xs0,2,hsPr0);
subplot(3,2,3); visualizeData(xs1,2,hs1);
subplot(3,2,4); visualizeData(xs1,2,hsPr1);
