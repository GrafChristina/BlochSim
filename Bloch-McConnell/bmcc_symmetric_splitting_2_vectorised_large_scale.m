function[M]=bmcc_symmetric_splitting_2_vectorised_large_scale(d)
%%% Symmetric Operator Splitting for 2 Pool classical CEST, large scale
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%

%A2 is constant in position and time
A2=zeros(6,6);
b=zeros(6,1);
b(3)=d.M0c(1)/d.T1(1)*d.relax(1);
b(6)=d.M0c(2)/d.T1(2)*d.relax(2);

A2(1,1)=-1/d.T2(1)*d.relax(1)-d.k(1);
A2(2,2)=-1/d.T2(1)*d.relax(1)-d.k(1);
A2(3,3)=-1/d.T1(1)*d.relax(1)-d.k(1);
A2(4,4)=-1/d.T2(2)*d.relax(2)-d.k(2);
A2(5,5)=-1/d.T2(2)*d.relax(2)-d.k(2);
A2(6,6)=-1/d.T1(2)*d.relax(2)-d.k(2);

A2(4,1)=d.k(1);
A2(5,2)=d.k(1);
A2(6,3)=d.k(1);
A2(1,4)=d.k(2);
A2(2,5)=d.k(2);
A2(3,6)=d.k(2);

if det(A2)==0
    disp('A2 not invertible! no relaxation and no coupling is assumend!');
    A2b=0*zeros(6,1);
else
    A2b=A2\b;
end


A_t1=[A2(1,1), A2(1,4); A2(4,1), A2(4,4)];
A_t3=[A2(3,3), A2(3,6); A2(6,3), A2(6,6)];

[V1, D1]=eig(A_t1);
V2=V1;
D2=D1;
[V3, D3]=eig(A_t3);
cases=zeros(1,3);
eps=10^-15;
va=[0;0];
vb=[0;0];
if abs(D1(1,1)-D1(2,2))<eps && abs(imag(D1(1,1)))<eps
    cases(1)=2;
    cases(2)=2;
    va=(A_t1-D1)\V1(:,1);
    va(isnan(va))=0;
    va(isinf(va))=0;
    if abs(A_t1(1,2))<eps && abs(A_t1(2,1))<eps
        cases(1)=3;
        cases(2)=3;
    end
elseif abs(imag(D1(1,1)))<eps
    cases(1)=3;
    cases(2)=3;
else
    cases(1)=1;
    cases(2)=1;
end

if abs(D3(1,1)-D3(2,2))<eps && abs(imag(D3(1,1)))<eps
    cases(3)=2;
    vb=(A_t3-D3)\V3(:,1);
    vb(isnan(vb))=0;
    vb(isinf(vb))=0;
    if abs(A_t3(1,2))<eps && abs(A_t3(2,1))<eps
        cases(3)=3;
        cases(3)=3;
    end
elseif abs(imag(D3(1,1))<eps)
    cases(3)=3;
else
    cases(3)=1;
end
Nx=d.Nx;
Nu=d.Nu;
M0=d.M0;
ud=d.u;

V=zeros(6,6);
V(1,1)=V1(1,1);
V(4,1)=V1(2,1);
V(1,4)=V1(1,2);
V(4,4)=V1(2,2);
V(2,2)=V2(1,1);
V(5,2)=V2(2,1);
V(2,5)=V2(1,2);
V(5,5)=V2(2,2);
V(3,3)=V3(1,1);
V(6,3)=V3(2,1);
V(3,6)=V3(1,2);
V(6,6)=V3(2,2);

D=[exp(D1(1,1)*d.dt)*V(:,1)*(cases(1)==3)+exp(real(D1(1,1))*d.dt)*(cos(imag(D1(1,1))*d.dt)*real(V(:,1))-sin(imag(D1(1,1))*d.dt)*imag(V(:,1)))*(cases(1)==1)+exp(D1(1,1)*d.dt)*V(:,1)*(cases(1)==2), ...
    exp(D1(2,2)*d.dt)*V(:,4)*(cases(1)==3)+exp(real(D1(1,1))*d.dt)*(sin(imag(D1(1,1))*d.dt)*real(V(:,1))+cos(imag(D1(1,1))*d.dt)*imag(V(:,1)))*(cases(1)==1)+(exp(D1(1,1)*d.dt)*d.dt*V(:,1)+[va(1);0;0;va(2);0;0])*(cases(1)==2), ...
    exp(D2(1,1)*d.dt)*V(:,2)*(cases(2)==3)+exp(real(D2(1,1))*d.dt)*(cos(imag(D2(1,1))*d.dt)*real(V(:,2))-sin(imag(D2(1,1))*d.dt)*imag(V(:,2)))*(cases(2)==1)+exp(D2(1,1)*d.dt)*V(:,2)*(cases(2)==2), ...
    exp(D2(2,2)*d.dt)*V(:,5)*(cases(2)==3)+exp(real(D2(1,1))*d.dt)*(sin(imag(D2(1,1))*d.dt)*real(V(:,2))+cos(imag(D2(1,1))*d.dt)*imag(V(:,2)))*(cases(2)==1)+(exp(D2(1,1)*d.dt)*d.dt*V(:,2)+[0;va(1);0;0;va(2);0])*(cases(2)==2), ...
    exp(D3(1,1)*d.dt)*V(:,3)*(cases(3)==3)+exp(real(D3(1,1))*d.dt)*(cos(imag(D3(1,1))*d.dt)*real(V(:,3))-sin(imag(D3(1,1))*d.dt)*imag(V(:,3)))*(cases(3)==1)+exp(D3(1,1).*d.dt)*V(:,3)*(cases(3)==2), ...
    exp(D3(2,2)*d.dt)*V(:,6)*(cases(3)==3)+exp(real(D3(1,1))*d.dt)*(sin(imag(D3(1,1))*d.dt)*real(V(:,3))+cos(imag(D3(1,1))*d.dt)*imag(V(:,3)))*(cases(3)==1)+(exp(D3(1,1)*d.dt)*d.dt*V(:,3)+[0;0;vb(1);0;0;vb(2)])*(cases(3)==2)];


detA=[(V1(1,1)*V1(2,2)-V1(2,1)*V1(1,2))*(cases(1)==3)+(real(V1(1,1))*imag(V1(2,2))-real(V1(2,1))*imag(V1(1,2)))*(cases(1)==1) + (V1(1,1)*va(2,1)-V1(2,1)*va(1,1))*(cases(1)==2), (V2(1,1)*V2(2,2)-V2(2,1)*V2(1,2))*(cases(2)==3)+(real(V2(1,1))*imag(V2(2,1))-real(V2(2,1))*imag(V2(1,1)))*(cases(2)==1)+(V2(1,1)*va(2,1)-V2(2,1)*va(1,1))*(cases(2)==2), (V3(1,1)*V3(2,2)-V3(2,1)*V3(1,2))*(cases(3)==3)+(real(V3(1,1))*imag(V3(2,1))-real(V3(2,1))*imag(V3(1,1)))*(cases(3)==1)+(V3(1,1)*vb(2,1)-V3(2,1)*vb(1,1))*(cases(3)==2)];

DC=zeros(6,6);
DC(1,1)=1/detA(1)*(V1(2,2)*(cases(1)==3)+imag(V1(2,2))*(cases(1)==1)+va(2,1)*(cases(1)==2));
DC(1,4)=1/detA(1)*(-V1(1,2)*(cases(1)==3)-imag(V1(1,2))*(cases(1)==1)-va(1,1)*(cases(1)==2));
DC(2,1)=1/detA(1)*(-V1(2,1)*(cases(1)==3)-real(V1(2,2))*(cases(1)==1)-V1(2,1)*(cases(1)==2));
DC(2,4)=1/detA(1)*(V1(1,1)*(cases(1)==3)+real(V1(1,2))*(cases(1)==1)+V1(1,1)*(cases(1)==2));
DC(3,2)=1/detA(2)*(V2(2,2)*(cases(2)==3)+imag(V1(2,2))*(cases(2)==1)+va(2,1)*(cases(2)==2));
DC(3,5)=1/detA(2)*(-V2(1,2)*(cases(2)==3)-imag(V1(1,2))*(cases(2)==1)-va(1,1)*(cases(2)==2));
DC(4,2)=1/detA(2)*(-V1(2,1)*(cases(2)==3)-real(V1(2,2))*(cases(2)==1)-V1(2,1)*(cases(2)==2));
DC(4,5)=1/detA(2)*(V1(1,1)*(cases(2)==3)+real(V1(1,2))*(cases(2)==1)+V1(1,1)*(cases(2)==2));
DC(5,3)=1/detA(3)*(V3(2,2)*(cases(3)==3)+imag(V3(2,2))*(cases(3)==1)+vb(2,1)*(cases(3)==2));
DC(5,6)=1/detA(3)*(-V3(1,2)*(cases(3)==3)-imag(V3(1,2))*(cases(3)==1)-vb(1,1)*(cases(3)==2));
DC(6,3)=1/detA(3)*(-V3(2,1)*(cases(3)==3)-real(V3(2,2))*(cases(3)==1)-V3(2,1)*(cases(3)==2));
DC(6,6)=1/detA(3)*(V3(1,1)*(cases(3)==3)+real(V3(1,2))*(cases(3)==1)+V3(1,1)*(cases(3)==2));
DC(isnan(DC))=0;
DC(isinf(DC))=0;


Dx=diag([cases(1)==4, cases(2)==4, cases(3)==4, cases(1)==4, cases(2)==4, cases(3)==4]);


d.u=ud;
M=zeros(6,d.Nx);

M(:,:)=M0;
Mt=M0;
d.u(d.u==0)=10^-14;
d.v(d.v==0)=10^-14;
for n=1:Nu-1 % time loop
    %%% Bloch simulation in magnetization domain
    gadt = d.gamma*d.dt/2;
    
    B = repmat(gadt*transpose(d.u(n)-1i*d.v(n))*d.B1c, Nx,1);
    wref=d.gamma*d.B0;
    w=wref*d.xZspec;
    K1=repmat((-w+wref*d.dw(1))*d.dt/2,1,1)';
    K2=repmat((-w+wref*d.dw(2))*d.dt/2,1,1)';
    phi1 = -(abs(B).^2+K1.^2).^(1/2);
    phi2 = -(abs(B).^2+K2.^2).^(1/2);
    
    cs1=cos(phi1);
    si1=sin(phi1);
    cs2=cos(phi2);
    si2=sin(phi2);
    
    n1 = real(B)./abs(phi1);
    n2 = imag(B)./abs(phi1);
    n3 = K1./abs(phi1);
    n4 = real(B)./abs(phi2);
    n5 = imag(B)./abs(phi2);
    n6 = K2./abs(phi2);
    
    
    Bd1 = n1.*n1.*(1-cs1)+cs1;
    Bd2 = n1.*n2.*(1-cs1)-n3.*si1;
    Bd3 = n1.*n3.*(1-cs1)+n2.*si1;
    Bd4 = n2.*n1.*(1-cs1)+n3.*si1;
    Bd5 = n2.*n2.*(1-cs1)+cs1;
    Bd6 = n2.*n3.*(1-cs1)-n1.*si1;
    Bd7 = n3.*n1.*(1-cs1)-n2.*si1;
    Bd8 = n3.*n2.*(1-cs1)+n1.*si1;
    Bd9 = n3.*n3.*(1-cs1)+cs1;
    
    Bd10 = n4.*n4.*(1-cs2)+cs2;
    Bd11 = n4.*n5.*(1-cs2)-n6.*si2;
    Bd12 = n4.*n6.*(1-cs2)+n5.*si2;
    Bd13 = n5.*n4.*(1-cs2)+n6.*si2;
    Bd14 = n5.*n5.*(1-cs2)+cs2;
    Bd15 = n5.*n6.*(1-cs2)-n4.*si2;
    Bd16 = n6.*n4.*(1-cs2)-n5.*si2;
    Bd17 = n6.*n5.*(1-cs2)+n4.*si2;
    Bd18 = n6.*n6.*(1-cs2)+cs2;
    
    %Rotation
    %             Mrot=Bd*Mn;
    Mrot(1,:)=Bd1'.*Mt(1,:)+Bd2'.*Mt(2,:)+Bd3'.*Mt(3,:);
    Mrot(2,:)=Bd4'.*Mt(1,:)+Bd5'.*Mt(2,:)+Bd6'.*Mt(3,:);
    Mrot(3,:)=Bd7'.*Mt(1,:)+Bd8'.*Mt(2,:)+Bd9'.*Mt(3,:);
    Mrot(4,:)=Bd10'.*Mt(4,:)+Bd11'.*Mt(5,:)+Bd12'.*Mt(6,:);
    Mrot(5,:)=Bd13'.*Mt(4,:)+Bd14'.*Mt(5,:)+Bd15'.*Mt(6,:);
    Mrot(6,:)=Bd16'.*Mt(4,:)+Bd17'.*Mt(5,:)+Bd18'.*Mt(6,:);
    
    Mrot_n=Mrot+A2b;
    
    C=DC*Mrot_n;
    Mt=D*C;
    
    Mt=Mt-A2b;
    
    
    Mtn=Dx*Mrot;
    
    Mt=Mt+Mtn;
    
    
    %Rotation
    Mn(1,:)=Bd1'.*Mt(1,:)+Bd2'.*Mt(2,:)+Bd3'.*Mt(3,:);
    Mn(2,:)=Bd4'.*Mt(1,:)+Bd5'.*Mt(2,:)+Bd6'.*Mt(3,:);
    Mn(3,:)=Bd7'.*Mt(1,:)+Bd8'.*Mt(2,:)+Bd9'.*Mt(3,:);
    Mn(4,:)=Bd10'.*Mt(4,:)+Bd11'.*Mt(5,:)+Bd12'.*Mt(6,:);
    Mn(5,:)=Bd13'.*Mt(4,:)+Bd14'.*Mt(5,:)+Bd15'.*Mt(6,:);
    Mn(6,:)=Bd16'.*Mt(4,:)+Bd17'.*Mt(5,:)+Bd18'.*Mt(6,:);
    
    
    Mt=Mn;
    M=Mn;
end
end