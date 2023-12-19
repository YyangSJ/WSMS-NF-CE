function FF_Dic = DFT_Dic(d,lambda,N,G)

FF_Dic=[];

for theta=-0.75+2/G:2/G:0.75
        r_s=[];G_s=[];
        for n=1:N 
            G_s(n)=1/sqrt(N)*exp(1i*2*pi/lambda*((n-1)*d*theta));
        end
        G_s=G_s.';
        FF_Dic=[FF_Dic G_s];
end
end