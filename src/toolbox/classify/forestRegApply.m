function [ys,ps, pd] = forestRegApply( data, forest, maxDepth, minCount, best )
% Apply learned forest classifier.
%
% USAGE
%  [ys,ps] = forestRegApply( data, forest, [maxDepth], [minCount], [best] )
%
% INPUTS
%  data     - [NxF] N length F feature vectors
%  forest   - learned forest classification model
%  maxDepth - [] maximum depth of tree
%  minCount - [] minimum number of data points to allow split
%  best     - [0] if true use single best prediction per tree
%
% OUTPUTS
%  ys       - [Nx1] predicted output labels
%  ps       - [NxH] predicted output label probabilities
%  pd       - [] predicted distribution
%
% EXAMPLE
%
% See also forestTrain
%
% Piotr's Image&Video Toolbox      Version 3.24
% Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see external/bsd.txt]
if(nargin<3 || isempty(maxDepth)), maxDepth=0; end
if(nargin<4 || isempty(minCount)), minCount=0; end
if(nargin<5 || isempty(best)), best=0; end
assert(isa(data,'single')); M=length(forest);
H=size(forest(1).distr,2); N=size(data,1);
mu = zeros(N,1);
variance = zeros(N,1);
if(best), ys=zeros(N,M); else ps=zeros(N,H); end
for i=1:M, tree=forest(i);
  if(maxDepth>0), tree.child(tree.depth>=maxDepth) = 0; end
  if(minCount>0), tree.child(tree.count<=minCount) = 0; end 
  ids = forestInds(data,tree.thrs,tree.fids,tree.child);
  if(best), ys(:,i)=tree.ys(ids); else ps=ps+tree.distr(ids,:); end
  mu = mu + tree.mean(ids);
  variance = variance + tree.var(ids);
end
if(best), ps=histc(ys',1:H)'; end; [~,ys]=max(ps,[],2); ps=ps/M;
mu = mu/M; variance = variance/M;
pd = [mu variance];
end

function [ids] = findRegInds(m, v, data)
%find index of classes via max of distributions
%m ... mean
%v ... variance
%data. data

end
