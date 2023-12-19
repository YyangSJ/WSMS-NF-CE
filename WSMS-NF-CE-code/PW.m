function g_s = PW(theta,d,lambda,N)

    for n=1:N
      %  r_s(m,n)=sqrt(r^2-2*r*((m-1)*D+(n-1)*d)*theta+(((m-1)*D+(n-1)*d))^2);  
        r_s(n)=-((n-1)*d)*theta;  
        % r_s(m,n)=r-((m-1)*D+(n-1)*d)*theta+(((m-1)*D))^2/2/r*(1-theta^2);  
        G_s(n)=1/sqrt(N)*exp(-1i*2*pi/lambda*(r_s(n))); 
    end

g_s=G_s(:);
end

