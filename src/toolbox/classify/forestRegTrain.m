function forest = forestRegTrain( data, ys, varargin )
% Train random forest classifier.
%
% Dimensions:
%  M - number trees
%  F - number features
%  N - number input vectors
%  H - number classes
%
% USAGE
%  forest = forestRegTrain( data, ys, [varargin] )
%
% INPUTS
%  data     - [NxF] N length F feature vectors
%  ys       - [Nx1] regressor output target -> m im paper
%  varargin - additional params (struct or name/value pairs)
%   .M          - [1] number of trees to train
%   .H          - [max(ys)] number of classes
%   .N1         - [5*N/M] number of data points for training each tree
%   .F1         - [sqrt(F)] number features to sample for each node split
%   .split      - ['gini'] options include 'gini', 'entropy' and 'twoing'
%   .minCount   - [1] minimum number of data points to allow split
%   .minChild   - [1] minimum number of data points allowed at child nodes
%   .maxDepth   - [64] maximum depth of tree
%   .dWts       - [] weights used for sampling and weighing each data point
%   .fWts       - [] weights used for sampling features
%
% OUTPUTS
%  forest   - learned forest model struct array w the following fields
%   .fids     - [Kx1] feature ids for each node
%   .thrs     - [Kx1] threshold corresponding to each fid
%   .child    - [Kx1] index of child for each node
%   .distr    - [KxH] prob distribution at each node
%   .ys       - [Kx1] or {Kx1} most likely label at each node
%   .count    - [Kx1] number of data points at each node
%   .depth    - [Kx1] depth of each node
%
% EXAMPLE
%  N=10000; H=5; d=2; [xs0,ys0,xs1,ys1]=demoGenData(N,N,H,d,1,1);
%  xs0=single(xs0); xs1=single(xs1);
%  pTrain={'maxDepth',50,'F1',2,'M',150,'minChild',5};
%  tic, forest=forestTrain(xs0,ys0,pTrain{:}); toc
%  ysPr0 = forestApply(xs0,forest);
%  ysPr1 = forestApply(xs1,forest);
%  e0=mean(ysPr0~=ys0); e1=mean(ysPr1~=ys1);
%  fprintf('errors trn=%f tst=%f\n',e0,e1); figure(1);
%  subplot(2,2,1); visualizeData(xs0,2,ys0);
%  subplot(2,2,2); visualizeData(xs0,2,ysPr0);
%  subplot(2,2,3); visualizeData(xs1,2,ys1);
%  subplot(2,2,4); visualizeData(xs1,2,ysPr1);
%
% See also forestApply, fernsClfTrain
%
% Piotr's Image&Video Toolbox      Version 3.24
% Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see external/bsd.txt]

% get additional parameters and fill in remaining parameters
dfs={ 'M',1, 'H',[], 'N1',[], 'F1',[], 'split','gini', 'minCount',1, ...
  'minChild',1, 'maxDepth',64, 'dWts',[], 'fWts',[] };
[M,H,N1,F1,splitStr,minCount,minChild,maxDepth,dWts,fWts] = ...
  getPrmDflt(varargin,dfs,1);
[~,F]=size(data); N = size(ys,1);  %do not make any error detection.....
minChild=max(1,minChild); minCount=max([1 minCount minChild]);
%if(isempty(H)), H=max(ys); end; assert(all(ys<=H)); %TODO: tidy up H
if(isempty(N1)), N1=round(5*N/M); end; N1=min(N,N1);
if(isempty(F1)), F1=round(sqrt(F)); end; F1=min(F,F1);
if(isempty(dWts)), dWts=ones(1,N,'single'); end; dWts=dWts/sum(dWts);
if(isempty(fWts)), fWts=ones(1,F,'single'); end; fWts=fWts/sum(fWts);
split=find(strcmpi(splitStr,{'gini','entropy','twoing','custom'}))-1;
if(isempty(split)), error('unknown splitting criteria: %s',splitStr); end

% make sure data has correct types
%if(~isa(data,'single')), data=single(data); end
%if(~isa(ys,'uint32')), ys=uint32(ys); end
if(~isa(fWts,'single')), fWts=single(fWts); end
if(~isa(dWts,'single')), dWts=single(dWts); end

% train M random trees on different subsets of data
prmTree = {H,F1,minCount,minChild,maxDepth,fWts,split};
for i=1:M
  if(N==N1), data1=data; ys1=ys; dWts1=dWts; else
    d=wswor(dWts,N1,4); data1=data(d,:); ys1=ys(d,:);
    dWts1=dWts(d); dWts1=dWts1/sum(dWts1);
  end
  tree = treeTrain(data1,ys1,dWts1,prmTree);
  if(mod(i,10)==0)
      i
  end
  if(i==1), forest=tree(ones(M,1)); else forest(i)=tree; end
end

end

function tree = treeTrain( data, ys, dWts, prmTree )
% Train single random tree.
[H,F1,minCount,minChild,maxDepth,fWts,split,]=deal(prmTree{:});
N=size(ys,1); %data size
K=2*N-1; %maximal number of nodes. E.g.: 2 nodes = 3

thrs=zeros(K,1,'single'); distr=zeros(K,H,'single');
%if(split ~= 3), distr=zeros(K,H,'single'); else distr = zeros(K, 2, 'single'); end
fids=zeros(K,1,'uint32'); child=fids; count=fids; depth=fids;
means=zeros(K,size(ys,2), 'double'); variances=zeros(K,size(ys,2), 'double');
d1=zeros(K, 1,'uint32'); d2=zeros(K,1,'uint32'); c1=zeros(K,1,'uint16'); c2=zeros(K,1,'uint16');
ysn=cell(K,1); dids=cell(K,1); dids{1}=uint32(1:N);
k=1; K=2; %k.. current node; K.. current number of nodes
while( k < K )
  dids1=dids{k}; dids{k}=[]; ys1=ys(dids1,:); n1=size(ys1,1); count(k)=n1;
  %distr(k,:)=histc(ys1,1:H)/n1; [~,ysn{k}]=max(distr(k,:));
  
  % if pure node or insufficient data don't train split
  if( n1<=minCount || depth(k)>maxDepth ), k=k+1; continue; end
  % train split and continue
  
  delta = randi(N,1,2);
  channel = randi(3,1,2)-1;  %channel: 0 - 2!!!
  
  fids1=wswor(fWts,F1,4);
  data1 = data{fids1}(delta(1), delta(2), channel(1), channel(2))'; %data1=data(dids1,fids1);
  data2 = data1(dids1);
  [~,order1]=sort(data2); order1=uint32(order1-1);
  [fid,thr,gain]=forestRegFindThr(data2,ys1,dWts(dids1),order1,split); %TODO: find splits, idee: mehrere zuf??llige werte, berechne objective function (entropie) f??r jeden, behalte den besten
  size(data2);
  fid=fids1(fid); left=data1(dids1)<thr; count0=nnz(left);
  
  if( gain>1e-10 && count0>=minChild && (n1-count0)>=minChild )
    child(k)=K; fids(k)=fid-1; thrs(k)=thr;
    dids{K}=dids1(left); dids{K+1}=dids1(~left);
    means(K,:)=mean(ys1(left,:)); means(K+1,:)=mean(ys1(~left,:));         %TODO: berechne gaussian der ??brig gebliebenen regression targets ys (m im paper)?
    variances(K,:)=var(double(ys1(left,:))); variances(K+1,:)=var(double(ys1(~left,:)));
    d1(k)=delta(1); d2(k)=delta(2); c1(k)=channel(1); c2(k)=channel(2);
    depth(K:K+1)=depth(k)+1; K=K+2;

  end; k=k+1;
end
% create output model struct
K=1:K-1; ysn=[ysn{K}]';
tree=struct('fids',fids(K),'thrs',thrs(K),'child',child(K),...
  'distr',distr(K,:),'ys',ysn,'count',count(K),'depth',depth(K),...
  'mean', means(K,:), 'var', variances(K,:),...
  'd1', d1(K), 'd2', d2(K), 'c1', c1(K), 'c2', c2(K));
%test if symoblic link is working as we expected
end

function ids = wswor( prob, N, trials )
% Fast weighted sample without replacement. Alternative to:
%  ids=datasample(1:length(prob),N,'weights',prob,'replace',false);
M=length(prob); assert(N<=M); if(N==M), ids=1:N; return; end
if(all(prob(1)==prob)), ids=randperm(M,N); return; end
cumprob=min([0 cumsum(prob)],1); assert(abs(cumprob(end)-1)<.01);
cumprob(end)=1; [~,ids]=histc(rand(N*trials,1),cumprob);
[s,ord]=sort(ids); K(ord)=[1; diff(s)]~=0; ids=ids(K);
if(length(ids)<N), ids=wswor(cumprob,N,trials*2); end
ids=ids(1:N)';
end
