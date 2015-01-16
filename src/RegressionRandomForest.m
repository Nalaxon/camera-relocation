%Autor: Thomas Pietsch, Bernhard Rapp
%Date:  23.10.2014
%Version 0.1
%Use Piotr's Matlab Toolbox [1] to train a regression random forest for camera
%relocation
%[1] http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

%do not forget load toolbox ;-)
%clear all;
clc; close all;

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
% N=1000; sig=.5; f=@(x) cos(x*pi*4+(x+1)).^2;
% xs0=rand(N,1); hs0=f(xs0)+randn(N,1)*sig;
% xs1=rand(N,1); hs1=f(xs1)+randn(N,1)*sig;
% xs0=single(xs0);
% xs1=single(xs1);


%% prepare data from files
color = imread('../data/heads/seq-01/frame-000000.color.png');
color = imresize(color, 2);
depth = imread('../data/heads/seq-01/frame-000000.depth.png');
depth = imresize(depth, 2);
pose = load('../data/heads/seq-01/frame-000000.pose.txt');
[a b] = size(depth);


D=@(p) int32(depth(p));
I=@(p) int32(color(p));
f_depth=@(d1,d2, c1, c2) D(d1:640*480+d1)-D(d2:640*480+d2); 
f_dargb=@(d1,d2,c1,c2) I(640*480*c1+d1:a*b*(c1+1)+d1) - I(640*480*c2+d2:a*b*(c2+1)+d2);
f_combined=@(d1,d2,c1,c2) f_dargb(d1,d2,c1,c2) + f_depth(d1,d2);

xs0 = repmat([{f_depth},{f_dargb},{f_combined}],a*b,1);

%for i=1:640
%    hs(1,640*(i-1)+1:640*i) = i;
%end
%hs(2,:) = repmat([1:480],1,640);

hs(1,:)=repmat([1:640],1,480);
for i=1:480
   hs(2,480*(i-1)+1:480*(i)) = i;
end
hs(3,:)=depth(1:640*480);
hs(4,:)=ones(1,640*480);

m = pose*hs;
m(1,:) = (m(1,:)./m(4,:));
m(2,:) = (m(2,:)./m(4,:));
m(3,:) = (m(3,:)./m(4,:));
m(4,:) = (m(4,:)./m(4,:));

hs0 = m(1:3,:)';

%hs0 = reshape(depth, [a*b 1]);

%% train
%train forest
pTrain={'maxDepth', 50, 'N1', a, 'F1', 1, 'M', 150, 'minChild', 1, 'split', 'custom'};

%train forst
forest=forestRegTrain(xs0, hs0, pTrain{:});


%% apply
%apply forst on first dataset
[hsPr0 ps0 pd0] = forestRegApply(xs0, forest);

%apply forst on second dataset
%[hsPr1 ps1 pd1] = forestRegApply(xs1, forest);

%some output they thought is usefull
e0=mean(hsPr0~=hs0); %e1=mean(hsPr1~=hs1);

sse=mean((hsPr0-hs0).^2);

fprintf('errors trn=%f sse=%f\n', e0, sse);
%fprintf('errors trn=%f tst=%f\n', e0, e1); figure(1);
%figure(1);
%subplot(2,2,1); visualizeData(xs0,2,hs0);
%subplot(2,2,2); visualizeData(xs0,2,hsPr0);
%subplot(3,2,3); visualizeData(xs1,2,hs1);
%subplot(3,2,4); visualizeData(xs1,2,hsPr1);
%t = 0:1:10000;

%subplot(2,2,3); plotDistribution(t, pd0);

%figure(2);
%figure(1); clf; hold on; plot(xs0,hs0,'.b'); plot(xs0,hsPr0,'.r'); hold off;
figure(1); clf; plot3(xs0(:,1), xs0(:,2), hs0, '.b'); hold on;
plot3(xs0(:,1), xs0(:,2), hsPr0, 'r.'); hold off;
title(sprintf('sse = %f', sse));

%subplot(3,2,6); plotDistribution(t, pd1, H);

