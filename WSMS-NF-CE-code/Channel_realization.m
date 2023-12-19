function [channel_SW,angle_grid,distance_grid,zz] = Channel_realization(d,D,lambda,L,M,N,G,angle_sample,distance_sample)
channel_SW=zeros(M*N,1);
 
z=1/sqrt(2)*(normrnd(0,1,L,1)+1i*normrnd(0,1,L,1));
angle_grid=randperm(length(angle_sample));
 
distance_grid=randperm(length(distance_sample));

for l=1:L
     theta(l)=angle_sample(angle_grid(l)); 
    r(l)=distance_sample(distance_grid(l)); 
    channel_SW=channel_SW+sqrt(M*N/L)*z(l)*SW(theta(l),r(l),d,D,lambda,M,N);
    zz(l)=sqrt(M*N/L)*z(l);
end

end

