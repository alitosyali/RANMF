function VV=GCSpectralClust2(A,Kmax,T)

N=length(A);
for n1=1:N;                                    % adding a few random edges improves performance
    for n2=n1+1:N; 
        if rand(1)<0.05; 
          A(n1,n2)=1; A(n2,n1)=1; 
        end; 
    end; 
end
D=zeros(N,N);
for n=1:N
	D(n,n)=sum(A(n,:));
end
AA=inv(D)*A;
[U,L]=eig(AA);
for K=1:Kmax
	[IDX,C]=kmeans(U(:,1:K),K,'EmptyAction','singleton','Start','uniform','Replicates',T); %,'Distance','cityblock');
	VV(:,K)=IDX;
end
