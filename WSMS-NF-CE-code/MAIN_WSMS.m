clear all;
close all;
rng(1)
f=100e9;
c=3e8;
lambda=c/f;
d=lambda/2;

M=8;N=24; %Number of RF chains and antennas
G= N; %angle grid resolution
D=lambda*8+N*d; %inter-subarray spacing
AP=D*M; % array aperture
R_NF=2*AP^2/lambda % near-field distance boundary
L=4;
L_est=L+1; % Estimated number of paths
P=8 % pilot overhead
angle_sample=-0.75+2/G:2/G:0.75;
distance_sample=5:5:2*R_NF; %users are located in far- and near-field grid
A=DFT_Dic(d,lambda,N,G);
B=NF_Dic_WSMS(d,D,lambda,N,M,angle_sample,distance_sample);
C=NF_Dic(D,lambda,M,G,angle_sample,distance_sample);
realization=5000;
SNR_sample=-15:5:15
for ii=1:realization
    ii
    for snr=1:length(SNR_sample)
        [h,theta,r,zz]=Channel_realization(d,D,lambda,L,M,N,G,angle_sample,distance_sample);
        W =normrnd(0,1,N,P)+1i*normrnd(0,1,N,P);
        W2 =normrnd(0,1,N,P)+1i*normrnd(0,1,N,P);
        [~,~,V_t]=svd(W./abs(W));
        [U_t,~,~]=svd(W2./abs(W2));
        W_opt=U_t*[eye(P);zeros(N-P,P)]*V_t'; 
        W_O=W_opt./(abs(W_opt));
%         W_G= W./abs(W);
          W_H=round(rand(N,P))*2-1;
%         W_I=[eye(P),zeros(P,N-P)]';
        Wt=W_O;
        Y_tilde=[]; Noi_tilde=[];W_tilde=[];
        for m=1:M
            for p=1:P
                Y(p,m)=Wt(:,p)'*h((m-1)*N+1:m*N)+Wt(:,p)'*sqrt(10^(-SNR_sample(snr)/10))*1/sqrt(2)*(normrnd(0,1,N,1)...
                    +1i*normrnd(0,1,N,1));
            end
        end
        for m=1:M
            W_tilde= blkdiag(W_tilde,Wt);
            Y_tilde=[Y_tilde;Y(:,m)];
        end
        %       for m=1:M
        %             Noi(:,m)=sqrt(10^(-SNR_sample(snr)/10))*1/sqrt(2)*(normrnd(0,1,P,1)+1i*normrnd(0,1,P,1));
        %             Y(:,m)=  Wt'*h((m-1)*N+1:m*N)+Noi(:,m);
        %             PHI(:,:,m)=Wt'*A;
        %       W_tilde= blkdiag(W_tilde,Wt);
        %             Noi_tilde=[Noi_tilde;Noi(:,m)];
        %    end
        %         Y_tilde=W_tilde'*h+Noi_tilde;
        

        %% PD-OMP 
        Gain_PD=cs_omp(vec(Y),W_tilde'*B,size(B,2),L_est);
        h_hat_PD=B*Gain_PD.';
        NMSE_PD(ii,snr)=norm(h-h_hat_PD(:))^2/norm(h)^2; 
        %% MAD-OMP
        h_hat_MAD=[];
        for m=1:M
            Gain_MAD =cs_omp(Y(:,m),Wt'*A,size(A,2),L_est);
            h_hat_MAD(:,m)=A*Gain_MAD.';
        end
        NMSE_MAD(ii,snr)=norm(h-h_hat_MAD(:))^2/norm(h)^2;
        
        
        %% TS-PAD-OMP 
        [Lambda,Gain] =cs_somp(Y,Wt'*A,size(A,2),L_est);
     
        theta_est=angle_sample(Lambda);
        Gain=Gain(Lambda,:);
        Channel_supp=[];Gain_new=[];
        
        for l=1:L_est
            %
            for rr=1:length(distance_sample)
                gain_r(rr)= conj(SW2(theta_est(l),distance_sample(rr),D,lambda,M))*Gain(l,:).';
                
            end
            
            [~,indmax]=max(abs(gain_r));
            r_est(l)=distance_sample(indmax);
            Channel_supp=[Channel_supp, SW(theta_est(l),r_est(l),d,D,lambda,M,N)];
            
        end
        G_est=W_tilde'*Channel_supp;
        h_hat_TS_PAD=Channel_supp*((G_est'*G_est)\G_est')*Y_tilde;
        NMSE_TS_PAD(ii,snr)=norm(h-h_hat_TS_PAD(:))^2/norm(h)^2; 
        Sup_true=[];WSup_true=[];
       %% 2D-PAD-OMP
        Sup=[];WSup=[]; 
        [h_hat_2D,theta_ind,r_ind,kappa]= x2D_PAD_OMP(Y,Wt,A,C,length(distance_sample),L_est);
  
        for l=1:L_est
            Sup=[Sup SW(angle_sample(theta_ind(l)),distance_sample(r_ind(l)),d,D,lambda,M,N)];
            WSup=[WSup W_tilde'*SW(angle_sample(theta_ind(l)),distance_sample(r_ind(l)),d,D,lambda,M,N)];
        end
        h_hat_2D=Sup*((WSup'*WSup)\WSup')*vec(Y);
        NMSE_2D_PAD(ii,snr)=norm(h-h_hat_2D)^2/norm(h)^2
        
        %% OLS lower bound
        for l=1:L
            Sup_true=[Sup_true SW(angle_sample(theta(l)),distance_sample(r(l)),d,D,lambda,M,N)];
            WSup_true=[WSup_true W_tilde'*SW(angle_sample(theta(l)),distance_sample(r(l)),d,D,lambda,M,N)];
        end
        h_hat_OLS=Sup_true*((WSup_true'*WSup_true)\WSup_true')*vec(Y);
        NMSE_OLS(ii,snr)=norm(h-h_hat_OLS)^2/norm(h)^2
    end
end

NMSE_PD_OUT=mean(NMSE_PD)
NMSE_MAD_OUT=mean(NMSE_MAD)
NMSE_TS_PAD_OUT=mean(NMSE_TS_PAD)
NMSE_2D_PAD2_OUT=mean(NMSE_2D_PAD)
NMSE_OLS_OUT=mean(NMSE_OLS)

co1= [0, 161, 241]/255;
co2=[29, 191, 151]/255
co3= [70, 158, 180]/255
co4=[253,185,106]/255
co5=[214,64,78]/255
figure
semilogy(SNR_sample, NMSE_PD_OUT, 'dk-', 'linewidth', 1.1, 'markerfacecolor', co1,'markersize', 7)
hold on
semilogy(SNR_sample, NMSE_MAD_OUT, '^k-', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 6.5)
hold on
semilogy(SNR_sample, NMSE_TS_PAD_OUT, 'sk-', 'linewidth', 1.1, 'markerfacecolor', co5,'markersize', 7.2)
hold on
semilogy(SNR_sample, NMSE_2D_PAD2_OUT, 'ok-', 'linewidth', 1.1, 'markerfacecolor', co4,'markersize', 6.5)
hold on
semilogy(SNR_sample, NMSE_OLS_OUT, '*k--', 'linewidth', 1.1, 'markerfacecolor', co2,'markersize', 7) 
hold on
grid on
 lgh=legend('PD-OMP','MAD-OMP','TS-PAD-OMP'...
     ,'2D-PAD-OMP','OLS');
 
set(lgh,'interpreter','latex', 'fontsize', 13); 
xlabel('SNR [dB]','interpreter','latex','fontsize',13)
ylabel('NMSE','interpreter','latex','fontsize',13)
 