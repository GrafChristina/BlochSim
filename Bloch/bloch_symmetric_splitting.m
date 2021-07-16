function[M]=bloch_symmetric_splitting(u,v,w,d)
%%% Symmetric splitting based solver for Bloch's equation with relaxation.
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%
xdis=d.xdis;                                                %spatial grid
Nx=d.Nx;                                                    %number of spatial points
Nu=size(u,1);                                               %number of time points
M0=d.M0;                                                    %initial magnetization 
dt=d.dt;                                                    %time step length

gadt = d.gamma*d.dt/2;                                      %time step multiplied with gamma
B1 = repmat(gadt*transpose(u-1i*v)*d.B1c, Nx,1);            %RF
K = gadt*xdis'*w'*d.G3;                                     %Gs

phi = -sqrt(abs(B1).^2+K.^2);                               %rotation angle
cs=cos(phi);
si=sin(phi);
n1 = real(B1)./abs(phi);                                    %rotation axis
n2 = imag(B1)./abs(phi);
n3 = K./abs(phi);
n1(isnan(n1))=1;
n2(isnan(n2))=0;
n3(isnan(n3))=0;
Bd1 = n1.*n1.*(1-cs)+cs;                                    %rotation matrix, 3x3
Bd2 = n1.*n2.*(1-cs)-n3.*si;
Bd3 = n1.*n3.*(1-cs)+n2.*si;
Bd4 = n2.*n1.*(1-cs)+n3.*si;
Bd5 = n2.*n2.*(1-cs)+cs;
Bd6 = n2.*n3.*(1-cs)-n1.*si;
Bd7 = n3.*n1.*(1-cs)-n2.*si;
Bd8 = n3.*n2.*(1-cs)+n1.*si;
Bd9 = n3.*n3.*(1-cs)+cs;


M=zeros(3, Nx,Nu);
M(:,:,1)=M0;
Mt=M0;
D=diag([exp(-1/d.T2*d.relax*dt), exp(-1/d.T2*d.relax*dt), exp(-1/d.T1*d.relax*dt)]); %relaxation
b=[0;0;d.M0c]-[0;0;d.M0c*exp(-1/d.T1*d.relax*dt)];          %right hand side


for n=1:Nu 
    %rotation
    Mrot(1,:)=Bd1(:,n)'.*Mt(1,:)+Bd2(:,n)'.*Mt(2,:)+Bd3(:,n)'.*Mt(3,:);
    Mrot(2,:)=Bd4(:,n)'.*Mt(1,:)+Bd5(:,n)'.*Mt(2,:)+Bd6(:,n)'.*Mt(3,:);
    Mrot(3,:)=Bd7(:,n)'.*Mt(1,:)+Bd8(:,n)'.*Mt(2,:)+Bd9(:,n)'.*Mt(3,:);
  
    %relaxation
    Mt=D*Mrot+b;
    
    %rotation
    Mrot(1,:)=Bd1(:,n)'.*Mt(1,:)+Bd2(:,n)'.*Mt(2,:)+Bd3(:,n)'.*Mt(3,:);
    Mrot(2,:)=Bd4(:,n)'.*Mt(1,:)+Bd5(:,n)'.*Mt(2,:)+Bd6(:,n)'.*Mt(3,:);
    Mrot(3,:)=Bd7(:,n)'.*Mt(1,:)+Bd8(:,n)'.*Mt(2,:)+Bd9(:,n)'.*Mt(3,:);

    
    Mt=Mrot;
    M(:,:,n+1)=Mrot;
end
end

