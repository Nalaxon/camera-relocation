function [ids] = forestRegInds(data, tree)
%UNTITLED2 finds the Indices of the forest
N = size(data, 1);

ids=zeros(N,1);

ids(:) = treeInds(data, ones(N,1), 1, tree, 1);


end

function [ids]=treeInds(data, dids, k, tree, depth)
if(tree.child(k)==0)
    ids(dids)=k;
    return;
end

depth
fid = tree.fids(k)+1;
data1 = data{fid}(tree.d1(k), tree.d2(k), tree.c1(k), tree.c2(k));

left = data1(dids) < tree.thrs(k);

ids(left) = treeInds(data, left, tree.child(k), tree, depth+1);

ids(~left) =treeInds(data, ~left, tree.child(k)+1, tree, depth+1);


end

