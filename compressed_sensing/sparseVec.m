function v=sparseVec(v,n)


tmp=sort(abs(v),'descend');
v(abs(v)<tmp(n))=0;