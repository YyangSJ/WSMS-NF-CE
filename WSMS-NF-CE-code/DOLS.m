function [Lambda,Gain] = DOLS(Y,PHI,L)
% Distributed orthogonal least squares for different sensing matrices
% sampling the common subspace.Written by Songjie Yang-9.13.2023.

[P,M]=size(Y);
I=1:size(PHI,2);
Lambda=[]; 
for k=1:L   
    for m=1:M
        for g=1:size(PHI,2) 
            if ismember(g,Lambda)
                continue;
            end
            PHI_G=PHI(:,[Lambda,g]);
            Res(g,m)=norm((eye(P)-PHI_G*((PHI_G'*PHI_G)\PHI_G'))*Y(:,m))^2/norm(Y(:,m))^2;
        end 
    end 
%     for m=1:M
%         Res_norm(:,m)=Res(:,m)/min(Res(:,m));
%     end
    Res_com=sum(Res,2);
    [~,indm]=min(Res_com);
    Lambda=[Lambda,indm];   
end


for m=1:M
    PHI_G=PHI(:,Lambda,m);
    Gain(:,m)=((PHI_G'*PHI_G)\PHI_G')*Y(:,m);
end
end
