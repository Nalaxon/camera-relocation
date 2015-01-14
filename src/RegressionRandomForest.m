%Autor: Thomas Pietsch, Bernhard Rapp
%Date:  23.10.2014
%Version 0.1
%Use Piotr's Matlab Toolbox [1] to train a regression random forest for camera
%relocation
%[1] http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

%do not forget load toolbox ;-)
clear; clc; close all;

toolboxRegCompile;


%define forst configuration
%N=10000; d=5; H=10;

%clc; clear all;
%[xs0, hs0, xs1, hs1] = demoGenData(N, N,H,d,1,1);
%M = csvread('../data/bike/day.csv');
%hs0 = M(:,16);
%xs0 = M(:,3:13);
%xs0=single(xs0); %xs1=single(xs1);


%% generate toy data
  N=1000; sig=.5; f=@(x) cos(x*pi*4)+(x+1).^2;
  xs0=rand(N,1); hs0=f(xs0)+randn(N,1)*sig;
  xs1=rand(N,1); hs1=f(xs1)+randn(N,1)*sig;

%train forest
pTrain={'maxDepth', 50, 'F1', 4 'M', 50, 'minChild', 1, 'split', 'custom'};

%train forst
forest=forestRegTrain(xs0, hs0, pTrain{:});

%apply forst on first dataset
[hsPr0 ps0 pd0] = forestRegApply(xs0, forest);

%apply forst on second dataset
%[hsPr1 ps1 pd1] = forestRegApply(xs1, forest);

%some output they thought is usefull
e0=mean(hsPr0~=hs0); %e1=mean(hsPr1~=hs1);

sse=mean((hsPr0-hs0).^2);

fprintf('errors trn=%f sse=%f\n', e0, sse); figure(1);
%fprintf('errors trn=%f tst=%f\n', e0, e1); figure(1);
figure(1);
subplot(2,2,1); visualizeData(xs0,2,hs0);
subplot(2,2,2); visualizeData(xs0,2,hsPr0);
%subplot(3,2,3); visualizeData(xs1,2,hs1);
%subplot(3,2,4); visualizeData(xs1,2,hsPr1);
t = 0:1:10000;

subplot(2,2,3); plotDistribution(t, pd0);

%subplot(3,2,6); plotDistribution(t, pd1, H);

