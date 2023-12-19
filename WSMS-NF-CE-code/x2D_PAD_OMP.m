function [h_est,theta_ind,r_ind,gain] = x2D_PAD_OMP(Y,W,A,C,r_number,L)

X=[];D_Gamma=[];
R=Y;
[~,A_number]=size(A);
[P,M]=size(Y);
for l=1:L
    for a=1:A_number
        TP(a,:)=abs(A(:,a)'*W*R*conj(C(:,(a-1)*r_number+1:a*r_number)))/norm(A(:,a)'*W)^2;
    end
    [x,y]=find(TP==max(max(TP)));
    theta_ind(l)=x(1);
    r_ind(l)=y(1);
    Dac=A(:,theta_ind(l))*C(:,(theta_ind(l)-1)*r_number+r_ind(l)).';
    D_Gamma=[D_Gamma, vec(Dac)];
    X=[X vec(W'*Dac)];
    kappa=((X'*X)\X')*vec(R);
    R=reshape(vec(R)-X*kappa,P,M);
end
gain=((X'*X)\X')*vec(Y);
h_est=D_Gamma*gain;
end

