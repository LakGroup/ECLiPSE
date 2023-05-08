function fd = LocalMinkowskiBouligandDim(A)

[n, r] = boxcount(A);

fd_list = -gradient(log(n))./gradient(log(r));
s_binary=reshape(find(diff([zeros(1,1) (abs(diff(fd_list))<0.1) zeros(1,1)])~=0),2,[]);
[lgtmax,jmax]=max(diff(s_binary));
istart=s_binary(1,jmax);
fd=mean(fd_list(istart:istart+lgtmax));

end