function Q=QFDistBased(V,A)
N=length(V);
K=max(V);
AV=zeros(N,N);
for i=1:K
	V1=find(V==i)';
	AV(V1,V1)=1;
end
Q=sum(sum(abs(A-AV)))/(N^2);
