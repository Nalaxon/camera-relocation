function [ids] = forestRegInds(data, tree)
%UNTITLED2 finds the Indices of the forest
N = size(data, 1);

%ids=zeros(N,1);

ids = treeInds(data, ones(N,1), 1, tree, 1);

end

function [ids]=treeInds(data, dids, k, tree, depth)
if(tree.child(k)==0)
    ids(dids)=k;
    return;
end

%depth


fid=tree.fids(k)+1;
data1=data{fid}(tree.d1(k), tree.d2(k), tree.c1(k), tree.c2(k));
left=find(data1(dids) < tree.thrs(k));
%size(left)
right=find(data1(dids) >= tree.thrs(k));
%size(right)

if(isempty(left))
    if(isempty(right)), return; else ids=zeros(max(right),1); end
else
    if(isempty(right)), ids=zeros(max(left),1); else ids=zeros(max(max(left),max(right)),1); end    
end

%size(ids)
if(~isempty(left))
t=treeInds(data, left, tree.child(k), tree, depth+1);
t(t==0)=[];
ids(left)=t;
end

if(~isempty(right))
t=treeInds(data, right, tree.child(k)+1, tree, depth+1);
t(t==0)=[];
ids(right)=t;
end
ids(ids==0)=[];
end

