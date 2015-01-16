function [ids] = forestRegInds(data, tree)
%UNTITLED2 finds the Indices of the forest
N = size(data, 1);
child = tree.child;
ids=zeros(N,1);
for i = 1:N
    k=1;
    while(child(k))
        fid = tree.fids(k);
        thr = tree.thrs(k);
        data1 = data{tree.fids(k)}(tree.d1(k), tree.d2(k), tree.c1(k), tree.c2(k));
        if(data1(i+fid*N) < thr)
            k=child(k)-1;
        else
            k=child(k);
        end
        ids(i)=k+1;
    end
end

end

