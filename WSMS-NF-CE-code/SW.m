function g_s = SW(theta,r,d,D,lambda,M,N)
for m=1:M
    for n=1:N
          r_s(m,n)=sqrt(r^2-2*r*((m-1)*D+(n-1)*d)*theta+(((m-1)*D+(n-1)*d))^2);  
       % r_s(m,n)=r-((m-1)*D+(n-1)*d)*theta+(((m-1)*D)+(n-1)*d)^2/2/r*(1-theta^2);  
     % r_s(m,n)=r-((m-1)*D+(n-1)*d)*theta+(((m-1)*D))^2/2/r*(1-theta^2);  
        G_s(n,m)=1/sqrt(M*N)*exp(-1i*2*pi/lambda*(r_s(m,n)-r)); 
    end
end
g_s=G_s(:);
end

