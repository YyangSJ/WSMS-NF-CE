function NF_Dic_WSMS = NF_Dic_WSMS(d,D,lambda,N,M,angle_sample,R_sample)

NF_Dic_WSMS=[];

for theta=angle_sample
    for r=R_sample
        r_s=[];G_s=[];
        for m=1:M
        for n=1:N
            r_s(n,m)=sqrt(r^2-2*r*((m-1)*D+(n-1)*d)*theta+(((m-1)*D+(n-1)*d))^2);  
            G_s(n,m)=1/sqrt(N*M)*exp(-1i*2*pi/lambda*(r_s(n,m)-r));
        end
        end
        G_s=G_s(:);
        NF_Dic_WSMS=[NF_Dic_WSMS G_s];
    end
end
end