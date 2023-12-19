function NF_Dic = NF_Dic(D,lambda,M,G,angle_sample,R_sample)

NF_Dic=[];

for theta=angle_sample
    for r=R_sample
        r_s=[];G_s=[];
        for m=1:M 
            r_s(m)=sqrt(r^2-2*r*((m-1)*D)*theta+(((m-1)*D))^2);  
            r_s(m)=r-((m-1)*D)*theta+(((m-1)*D))^2/2/r*(1-theta^2);  
            G_s(m)=1/sqrt(M)*exp(-1i*2*pi/lambda*(r_s(m)-r));
        end
        G_s=G_s(:);
        NF_Dic=[NF_Dic G_s];
    end
end
end