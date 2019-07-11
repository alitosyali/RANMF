function W=PermMat(N)

W=zeros(N,N);
q=randperm(N);
for n=1:N; 
	W(q(n),n)=1; 
end
