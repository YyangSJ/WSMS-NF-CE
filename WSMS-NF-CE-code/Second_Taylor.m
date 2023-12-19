function g_f = Second_Taylor(theta,phi,r,d,lambda,M,N)
for m=1:M
    for n=1:N 
        r_f(m,n)=r-m*d*theta-n*d*phi+(n*d)^2/2/r*(1-phi^2)+(m*d)^2/2/r*(1-theta^2)-m*n*d^2*theta*phi/r; 
        G_f(n,m)=1/sqrt(M*N)*exp(-1i*2*pi/lambda*(r_f(m,n)-r)); 
    end
end
g_f=G_f(:);
end

