function[M]=bmcc_symmetric_splitting_3_vectorised_large_scale(d)
%%% Symmetric Operator for 3 Pool classical CEST, large scale variant
% See
% C. Graf, A. Rund, C.S. Aigner, R. Stollberger,
% Accuracy and Performance Analysis for Bloch and Bloch-McConnell
% simulation methods
% Journal of Magnetic Resonance 329(3):107011
% doi: 10.1016/j.jmr.2021.107011
%%%

%%% Relaxation and transfer rates
A2=zeros(9,9);
b=zeros(9,1);
b(3)=d.M0c(1)/d.T1(1)*d.relax(1);
b(6)=d.M0c(2)/d.T1(2)*d.relax(2);
b(9)=d.M0c(3)/d.T1(3)*d.relax(3);

A2(1,1)=-1/d.T2(1)*d.relax(1)-d.k(1)-d.k(2);
A2(2,2)=-1/d.T2(1)*d.relax(1)-d.k(1)-d.k(2);
A2(3,3)=-1/d.T1(1)*d.relax(1)-d.k(1)-d.k(2);
A2(4,4)=-1/d.T2(2)*d.relax(2)-d.k(3)-d.k(4);
A2(5,5)=-1/d.T2(2)*d.relax(2)-d.k(3)-d.k(4);
A2(6,6)=-1/d.T1(2)*d.relax(2)-d.k(3)-d.k(4);
A2(7,7)=-1/d.T2(3)*d.relax(3)-d.k(5)-d.k(6);
A2(8,8)=-1/d.T2(3)*d.relax(3)-d.k(5)-d.k(6);
A2(9,9)=-1/d.T1(3)*d.relax(3)-d.k(5)-d.k(6);

A2(4,1)=d.k(1);
A2(5,2)=d.k(1);
A2(6,3)=d.k(1);
A2(1,4)=d.k(3);
A2(2,5)=d.k(3);
A2(3,6)=d.k(3);
A2(7,1)=d.k(2);
A2(8,2)=d.k(2);
A2(9,3)=d.k(2);
A2(1,7)=d.k(5);
A2(2,8)=d.k(5);
A2(3,9)=d.k(5);
A2(7,4)=d.k(4);
A2(8,5)=d.k(4);
A2(9,6)=d.k(4);
A2(4,7)=d.k(6);
A2(5,8)=d.k(6);
A2(6,9)=d.k(6);
A2b=A2\b;
A2b(isnan(A2b))=0;

A21=A2(1:3:end,1:3:end);
A23=A2(3:3:end,3:3:end);

[V1, D1]=eig(A21);
V2=V1;
D2=D1;
[V3, D3]=eig(A23);

cases=zeros(3,1);
eps=10^-15;
vc_1=[0;0;0];
vc_3=[0;0;0];

if norm(imag(D1))<eps
    if abs(D1(1,1))<eps && abs(D1(2,2))<eps && abs(D1(3,3))<eps
        cases(1)=2;
        cases(2)=2;
        lambda1=D1(1,1);
        lambda2=D1(2,2);
        lambda3=D1(3,3);
        Sk1=V1;
        lambda4=lambda1;
        lambda5=lambda2;
        lambda6=lambda3;
        Sk2=V2;
    elseif ((abs(D1(1,1)-D1(2,2))<eps && abs(D1(1,1)-D1(3,3))>eps) || (abs(D1(2,2)-D1(3,3))<eps && abs(D1(1,1)-D1(2,2))>eps)) && norm(d.k-[0 0 0 0 0 0])~=0
        cases(1)=3;
        cases(2)=3;
        vc_1=(A21-D1)\V1(:,1);
        vc_1(isnan(vc_1))=0;
        if abs(D1(1,1)-D1(2,2))<eps
            lambda1=D1(1,1);
            lambda2=D1(2,2);
            lambda3=D1(3,3);
            Sk1=V1;
            lambda4=lambda1;
            lambda5=lambda2;
            lambda6=lambda3;
            Sk2=V2;
        else
            lambda2=D1(1,1);
            lambda3=D1(2,2);
            lambda1=D1(3,3);
            Sk1=[V1(:,2) V1(:,3) V1(:,1)];
            lambda5=D2(1,1);
            lambda6=D2(2,2);
            lambda4=D2(3,3);
            Sk2=[V2(:,2) V2(:,3) V2(:,1)];
        end
    else
        cases(1)=4;
        cases(2)=4;
        lambda1=D1(1,1);
        lambda2=D1(2,2);
        lambda3=D1(3,3);
        Sk1=V1;
        lambda4=lambda1;
        lambda5=lambda2;
        lambda6=lambda3;
        Sk2=V2;
    end
else
    %%% a real eigenvalue and a complex conjugate pair, ordered s.t. the
    %%% real is the first one
    cases(1)=1;
    cases(2)=1;
    lambda1=D1(1,1);
    lambda2=D1(2,2);
    lambda3=D1(3,3);
    lambda4=lambda1;
    lambda5=lambda2;
    lambda6=lambda3;
    Sk1=V1;
    Sk2=V2;
    if abs(imag(lambda1))>eps
        lambda1=D1(3,3);
        lambda2=D1(1,1);
        lambda3=D1(2,2);
        lambda4=lambda1;
        lambda5=lambda2;
        lambda6=lambda3;
        Sk1=[V1(:,3) V1(:,1) V1(:,2)];
        Sk2=[V2(:,3) V2(:,1) V2(:,2)];
    end
end

if norm(imag(D3))<eps
    if abs(D3(1,1))<eps && abs(D3(2,2))<eps && abs(D3(3,3))<eps
        cases(3)=2;
        lambda7=D3(1,1);
        lambda8=D3(2,2);
        lambda9=D3(3,3);
        Sk3=V3;
    elseif ((abs(D3(1,1)-D3(2,2))<eps && abs(D3(1,1)-D3(3,3))>eps) || (abs(D3(2,2)-D3(3,3))<eps && abs(D3(1,1)-D3(2,2))>eps))&& norm(d.k-[0 0 0 0 0 0])~=0
        cases(3)=3;
        vc_3=(A23-D3)\V3(:,1);
        vc_3(isnan(vc_3))=0;
        if abs(D3(1,1)-D3(2,2))<eps
            lambda7=D3(1,1);
            lambda8=D3(2,2);
            lambda9=D3(3,3);
            Sk3=V3;
        else
            lambda7=D3(3,3);
            lambda8=D3(1,1);
            lambda9=D3(2,2);
            Sk3=[V3(:,2) V3(:,3) V3(:,1)];
        end
    else
        cases(3)=4;
        lambda7=D3(1,1);
        lambda8=D3(2,2);
        lambda9=D3(3,3);
        Sk3=V3;
    end
else
    %%% a real eigenvalue and a complex conjugate pair, ordered s.t. the
    %%% real is the first one
    cases(3)=1;
    lambda7=D3(1,1);
    lambda8=D3(2,2);
    lambda9=D3(3,3);
    Sk3=V3;
    if abs(imag(lambda7))>eps
        lambda7=D3(3,3);
        lambda8=D3(1,1);
        lambda9=D3(2,2);
        Sk3=[V3(:,3) V3(:,1) V3(:,2)];
    end
end



%%% eigenvectors of the whole matrix A2
V=zeros(9,9);
V(1,1)=Sk1(1,1);
V(4,1)=Sk1(2,1);
V(7,1)=Sk1(3,1);
V(1,4)=Sk1(1,2);
V(4,4)=Sk1(2,2);
V(7,4)=Sk1(3,2);
V(1,7)=Sk1(1,3);
V(4,7)=Sk1(2,3);
V(7,7)=Sk1(3,3);
V(2,2)=Sk2(1,1);
V(5,2)=Sk2(2,1);
V(8,2)=Sk2(3,1);
V(2,5)=Sk2(1,2);
V(5,5)=Sk2(2,2);
V(8,5)=Sk2(3,2);
V(2,8)=Sk2(1,3);
V(5,8)=Sk2(2,3);
V(8,8)=Sk2(3,3);
V(3,3)=Sk3(1,1);
V(6,3)=Sk3(2,1);
V(9,3)=Sk3(3,1);
V(3,6)=Sk3(1,2);
V(6,6)=Sk3(2,2);
V(9,6)=Sk3(3,2);
V(3,6)=Sk3(1,2);
V(6,6)=Sk3(2,2);
V(9,6)=Sk3(3,2);
V(3,9)=Sk3(1,3);
V(6,9)=Sk3(2,3);
V(9,9)=Sk3(3,3);

%%% generalized eigenvectors if no basis consisting of eigenvectors is
%%% available
vc_1=[vc_1(1);0;0;vc_1(2);0;0;vc_1(3);0;0];
vc_2=[0;vc_1(1);0;0;vc_1(2);0;0;vc_1(3);0];
vc_3=[0;0;vc_3(1);0;0;vc_3(2);0;0;vc_3(3)];


%%% solution basis
D=[V(:,1).*exp(lambda1*d.dt).*(cases(1)==1)+V(:,1)*exp(lambda1*d.dt).*(cases(1)==3)+V(:,1)*exp(lambda1*d.dt).*(cases(1)==4), ...
    V(:,2).*exp(lambda4*d.dt).*(cases(2)==1)+V(:,2)*exp(lambda4*d.dt).*(cases(2)==3)+V(:,2)*exp(lambda4*d.dt).*(cases(2)==4), ...
    V(:,3).*exp(lambda7*d.dt).*(cases(3)==1)+V(:,3)*exp(lambda7*d.dt).*(cases(3)==3)+V(:,3)*exp(lambda7*d.dt).*(cases(3)==4), ...
    exp(real(lambda2)*d.dt)*(cos(imag(lambda2)*d.dt)*real(V(:,4))-sin(imag(lambda2)*d.dt)*imag(V(:,4))).*(cases(1)==1)+exp(lambda2*d.dt)*(vc_1+d.dt*V(:,4)).*(cases(1)==3)+V(:,4)*exp(lambda2*d.dt).*(cases(1)==4), ...
    exp(real(lambda5)*d.dt)*(cos(imag(lambda5)*d.dt)*real(V(:,5))-sin(imag(lambda5)*d.dt)*imag(V(:,5))).*(cases(2)==1)+exp(lambda5*d.dt)*(vc_1+d.dt*V(:,5)).*(cases(2)==3)+V(:,5)*exp(lambda5*d.dt).*(cases(2)==4), ...
    exp(real(lambda8)*d.dt)*(cos(imag(lambda8)*d.dt)*real(V(:,6))-sin(imag(lambda8)*d.dt)*imag(V(:,6))).*(cases(3)==1)+exp(lambda8*d.dt)*(vc_3+d.dt*V(:,6)).*(cases(3)==3)+V(:,6)*exp(lambda8*d.dt).*(cases(3)==4), ...
    exp(real(lambda2)*d.dt)*(sin(imag(lambda2)*d.dt)*real(V(:,4))+cos(imag(lambda2)*d.dt)*imag(V(:,4))).*(cases(1)==1)+V(:,7)*exp(lambda3*d.dt).*(cases(1)==3)+V(:,7)*exp(lambda3*d.dt).*(cases(1)==4), ...
    exp(real(lambda5)*d.dt)*(sin(imag(lambda5)*d.dt)*real(V(:,5))+cos(imag(lambda5)*d.dt)*imag(V(:,5))).*(cases(2)==1)+V(:,8)*exp(lambda6*d.dt).*(cases(2)==3)+V(:,8)*exp(lambda6*d.dt).*(cases(2)==4), ...
    exp(real(lambda8)*d.dt)*(sin(imag(lambda8)*d.dt)*real(V(:,6))+cos(imag(lambda8)*d.dt)*imag(V(:,6))).*(cases(3)==1)+V(:,9)*exp(lambda9*d.dt).*(cases(3)==3)+V(:,9)*exp(lambda9*d.dt).*(cases(3)==4)];

vc_1=[vc_1(1);vc_1(4);vc_1(7)];
vc_2=vc_1;
vc_3=[vc_3(3);vc_3(6);vc_3(9)];
detA=[det([Sk1(:,1) real(Sk1(:,2)) imag(Sk1(:,2))]).*(cases(1)==1)+det([Sk1(:,1),Sk1(:,3),vc_1]).*(cases(1)==3)+det([Sk1(:,1) Sk1(:,2) Sk1(:,3)]).*(cases(1)==4), det([Sk2(:,1) real(Sk2(:,2)) imag(Sk2(:,2))]).*(cases(2)==1)+det([Sk2(:,1),Sk2(:,3),vc_1]).*(cases(2)==3)+det([Sk2(:,1) Sk2(:,2) Sk2(:,3)]).*(cases(2)==4),det([Sk3(:,1) real(Sk3(:,2)) imag(Sk3(:,2))]).*(cases(3)==1)+det([Sk3(:,1),Sk3(:,3),vc_3]).*(cases(3)==3)+det([Sk3(:,1) Sk3(:,2) Sk3(:,3)]).*(cases(3)==4)];


%%% auxiliary matrix for the calculation of the constants
DC=zeros(9,9);
DC(1,1)=1/detA(1)*((real(Sk1(2,2))*imag(Sk1(3,2))-real(Sk1(3,2))*imag(Sk1(2,2))).*(cases(1)==1)+(Sk1(2,3)*vc_1(3)-Sk1(3,3)*vc_1(2)).*(cases(1)==3)+(Sk1(2,2)*Sk1(3,3)-Sk1(3,2)*Sk1(2,3)).*(cases(1)==4));
DC(1,4)=1/detA(1)*((real(Sk1(3,2))*imag(Sk1(1,2))-real(Sk1(1,2))*imag(Sk1(3,2))).*(cases(1)==1)+(Sk1(3,3)*vc_1(1)-Sk1(1,3)*vc_1(3)).*(cases(1)==3)+(Sk1(3,2)*Sk1(1,3)-Sk1(1,2)*Sk1(3,3)).*(cases(1)==4));
DC(1,7)=1/detA(1)*((real(Sk1(1,2))*imag(Sk1(2,2))-real(Sk1(2,2))*imag(Sk1(1,2))).*(cases(1)==1)+(Sk1(1,3)*vc_1(2)-Sk1(2,3)*vc_1(1)).*(cases(1)==3)+(Sk1(1,2)*Sk1(2,3)-Sk1(2,2)*Sk1(1,3)).*(cases(1)==4));
DC(4,1)=1/detA(1)*((Sk1(3,1)*imag(Sk1(2,2))-Sk1(2,1)*imag(Sk1(3,2))).*(cases(1)==1)+(Sk1(3,1)*vc_1(2)-Sk1(2,1)*vc_1(3)).*(cases(1)==3)+(Sk1(3,1)*Sk1(2,3)-Sk1(2,1)*Sk1(3,3)).*(cases(1)==4));
DC(4,4)=1/detA(1)*((Sk1(1,1)*imag(Sk1(3,2))-Sk1(3,1)*imag(Sk1(1,2))).*(cases(1)==1)+(Sk1(1,1)*vc_1(3)-Sk1(3,1)*vc_1(1)).*(cases(1)==3)+(Sk1(1,1)*Sk1(3,3)-Sk1(3,1)*Sk1(1,3)).*(cases(1)==4));
DC(4,7)=1/detA(1)*((Sk1(2,1)*imag(Sk1(1,2))-Sk1(1,1)*imag(Sk1(2,2))).*(cases(1)==1)+(Sk1(2,1)*vc_1(1)-Sk1(1,1)*vc_1(2)).*(cases(1)==3)+(Sk1(2,1)*Sk1(1,3)-Sk1(1,1)*Sk1(2,3)).*(cases(1)==4));
DC(7,1)=1/detA(1)*((Sk1(2,1)*real(Sk1(3,2))-Sk1(3,1)*real(Sk1(2,2))).*(cases(1)==1)+(Sk1(2,1)*Sk1(3,3)-Sk1(3,1)*Sk1(2,3)).*(cases(1)==3)+(Sk1(2,1)*Sk1(3,2)-Sk1(3,1)*Sk1(2,2)).*(cases(1)==4));
DC(7,4)=1/detA(1)*((Sk1(3,1)*real(Sk1(1,2))-Sk1(1,1)*real(Sk1(3,2))).*(cases(1)==1)+(Sk1(3,1)*Sk1(1,3)-Sk1(1,1)*Sk1(3,3)).*(cases(1)==3)+(Sk1(3,1)*Sk1(1,2)-Sk1(1,1)*Sk1(3,2)).*(cases(1)==4));
DC(7,7)=1/detA(1)*((Sk1(1,1)*real(Sk1(2,2))-Sk1(2,1)*real(Sk1(1,2))).*(cases(1)==1)+(Sk1(1,1)*Sk1(2,3)-Sk1(2,1)*Sk1(1,3)).*(cases(1)==3)+(Sk1(1,1)*Sk1(2,2)-Sk1(2,1)*Sk1(1,2)).*(cases(1)==4));

DC(2,2)=1/detA(2)*((real(Sk2(2,2))*imag(Sk2(3,2))-real(Sk2(3,2))*imag(Sk2(2,2))).*(cases(2)==1)+(Sk2(2,3)*vc_1(3)-Sk2(3,3)*vc_1(2)).*(cases(2)==3)+(Sk2(2,2)*Sk2(3,3)-Sk2(3,2)*Sk2(2,3)).*(cases(2)==4));
DC(2,5)=1/detA(2)*((real(Sk2(3,2))*imag(Sk2(1,2))-real(Sk2(1,2))*imag(Sk2(3,2))).*(cases(2)==1)+(Sk2(3,3)*vc_1(1)-Sk2(1,3)*vc_1(3)).*(cases(2)==3)+(Sk2(3,2)*Sk2(1,3)-Sk2(1,2)*Sk2(3,3)).*(cases(2)==4));
DC(2,8)=1/detA(2)*((real(Sk2(1,2))*imag(Sk2(2,2))-real(Sk2(2,2))*imag(Sk2(1,2))).*(cases(2)==1)+(Sk2(1,3)*vc_1(2)-Sk2(2,3)*vc_1(1)).*(cases(2)==3)+(Sk2(1,2)*Sk2(2,3)-Sk2(2,2)*Sk2(1,3)).*(cases(2)==4));
DC(5,2)=1/detA(2)*((Sk2(3,1)*imag(Sk2(2,2))-Sk2(2,1)*imag(Sk2(3,2))).*(cases(2)==1)+(Sk2(3,1)*vc_1(2)-Sk2(2,1)*vc_1(3)).*(cases(2)==3)+(Sk2(3,1)*Sk2(2,3)-Sk2(2,1)*Sk2(3,3)).*(cases(2)==4));
DC(5,5)=1/detA(2)*((Sk2(1,1)*imag(Sk2(3,2))-Sk2(3,1)*imag(Sk2(1,2))).*(cases(2)==1)+(Sk2(1,1)*vc_1(3)-Sk2(3,1)*vc_1(1)).*(cases(2)==3)+(Sk2(1,1)*Sk2(3,3)-Sk2(3,1)*Sk2(1,3)).*(cases(2)==4));
DC(5,8)=1/detA(2)*((Sk2(2,1)*imag(Sk2(1,2))-Sk2(1,1)*imag(Sk2(2,2))).*(cases(2)==1)+(Sk2(2,1)*vc_1(1)-Sk2(1,1)*vc_1(2)).*(cases(2)==3)+(Sk2(2,1)*Sk2(1,3)-Sk2(1,1)*Sk2(2,3)).*(cases(2)==4));
DC(8,2)=1/detA(2)*((Sk2(2,1)*real(Sk2(3,2))-Sk2(3,1)*real(Sk2(2,2))).*(cases(2)==1)+(Sk2(2,1)*Sk2(3,3)-Sk2(3,1)*Sk2(2,3)).*(cases(2)==3)+(Sk2(2,1)*Sk2(3,2)-Sk2(3,1)*Sk2(2,2)).*(cases(2)==4));
DC(8,5)=1/detA(2)*((Sk2(3,1)*real(Sk2(1,2))-Sk2(1,1)*real(Sk2(3,2))).*(cases(2)==1)+(Sk2(3,1)*Sk2(1,3)-Sk2(1,1)*Sk2(3,3)).*(cases(2)==3)+(Sk2(3,1)*Sk2(1,2)-Sk2(1,1)*Sk2(3,2)).*(cases(2)==4));
DC(8,8)=1/detA(2)*((Sk2(1,1)*real(Sk2(2,2))-Sk2(2,1)*real(Sk2(1,2))).*(cases(2)==1)+(Sk2(1,1)*Sk2(2,3)-Sk2(2,1)*Sk2(1,3)).*(cases(2)==3)+(Sk2(1,1)*Sk2(2,2)-Sk2(2,1)*Sk2(1,2)).*(cases(2)==4));

DC(3,3)=1/detA(3)*((real(Sk3(2,2))*imag(Sk3(3,2))-real(Sk3(3,2))*imag(Sk3(2,2))).*(cases(3)==1)+(Sk3(2,3)*vc_3(3)-Sk3(3,3)*vc_3(2)).*(cases(3)==3)+(Sk3(2,2)*Sk3(3,3)-Sk3(3,2)*Sk3(2,3)).*(cases(3)==4));
DC(3,6)=1/detA(3)*((real(Sk3(3,2))*imag(Sk3(1,2))-real(Sk3(1,2))*imag(Sk3(3,2))).*(cases(3)==1)+(Sk3(3,3)*vc_3(1)-Sk3(1,3)*vc_3(3)).*(cases(3)==3)+(Sk3(3,2)*Sk3(1,3)-Sk3(1,2)*Sk3(3,3)).*(cases(3)==4));
DC(3,9)=1/detA(3)*((real(Sk3(1,2))*imag(Sk3(2,2))-real(Sk3(2,2))*imag(Sk3(1,2))).*(cases(3)==1)+(Sk3(1,3)*vc_3(2)-Sk3(2,3)*vc_3(1)).*(cases(3)==3)+(Sk3(1,2)*Sk3(2,3)-Sk3(2,2)*Sk3(1,3)).*(cases(3)==4));
DC(6,3)=1/detA(3)*((Sk3(3,1)*imag(Sk3(2,2))-Sk3(2,1)*imag(Sk3(3,2))).*(cases(3)==1)+(Sk3(3,1)*vc_3(2)-Sk3(2,1)*vc_3(3)).*(cases(3)==3)+(Sk3(3,1)*Sk3(2,3)-Sk3(2,1)*Sk3(3,3)).*(cases(3)==4));
DC(6,6)=1/detA(3)*((Sk3(1,1)*imag(Sk3(3,2))-Sk3(3,1)*imag(Sk3(1,2))).*(cases(3)==1)+(Sk3(1,1)*vc_3(3)-Sk3(3,1)*vc_3(1)).*(cases(3)==3)+(Sk3(1,1)*Sk3(3,3)-Sk3(3,1)*Sk3(1,3)).*(cases(3)==4));
DC(6,9)=1/detA(3)*((Sk3(2,1)*imag(Sk3(1,2))-Sk3(1,1)*imag(Sk3(2,2))).*(cases(3)==1)+(Sk3(2,1)*vc_3(1)-Sk3(1,1)*vc_3(2)).*(cases(3)==3)+(Sk3(2,1)*Sk3(1,3)-Sk3(1,1)*Sk3(2,3)).*(cases(3)==4));
DC(9,3)=1/detA(3)*((Sk3(2,1)*real(Sk3(3,2))-Sk3(3,1)*real(Sk3(2,2))).*(cases(3)==1)+(Sk3(2,1)*Sk3(3,3)-Sk3(3,1)*Sk3(2,3)).*(cases(3)==3)+(Sk3(2,1)*Sk3(3,2)-Sk3(3,1)*Sk3(2,2)).*(cases(3)==4));
DC(9,6)=1/detA(3)*((Sk3(3,1)*real(Sk3(1,2))-Sk3(1,1)*real(Sk3(3,2))).*(cases(3)==1)+(Sk3(3,1)*Sk3(1,3)-Sk3(1,1)*Sk3(3,3)).*(cases(3)==3)+(Sk3(3,1)*Sk3(1,2)-Sk3(1,1)*Sk3(3,2)).*(cases(3)==4));
DC(9,9)=1/detA(3)*((Sk3(1,1)*real(Sk3(2,2))-Sk3(2,1)*real(Sk3(1,2))).*(cases(3)==1)+(Sk3(1,1)*Sk3(2,3)-Sk3(2,1)*Sk3(1,3)).*(cases(3)==3)+(Sk3(1,1)*Sk3(2,2)-Sk3(2,1)*Sk3(1,2)).*(cases(3)==4));

DC(isnan(DC))=0;

Dx=zeros(9,9);
Dx(1:1:3,1:1:3)=eye(3).*(cases(1)==2);
Dx(4:1:6,4:1:6)=eye(3).*(cases(2)==2);
Dx(7:1:9,7:1:9)=eye(3).*(cases(3)==2);




M=zeros(9,d.Nx);
M(:,:)=d.M0;
Mt=d.M0;
d.u(d.u==0)=10^-14;
d.v(d.v==0)=10^-14;
for n=1:d.Nu-1 % time loop
    
    %%% Bloch simulation in magnetization domain
    gadt = d.gamma*d.dt/2;
    d.u(d.u==0)=10^-14;
    d.v(d.v==0)=10^-14;
    B = repmat(gadt*transpose(d.u(n)-1i*d.v(n))*d.B1c, d.Nx,1);
    wref=d.gamma*d.B0;
    w=wref*d.xZspec;
    K1=repmat((-w+wref*d.dw(1))*d.dt/2,1,1)';
    K2=repmat((-w+wref*d.dw(2))*d.dt/2,1,1)';
    K3=repmat((-w+wref*d.dw(3))*d.dt/2,1,1)';
    
    phi1 = -sqrt(abs(B).^2+K1.^2);
    phi2 = -sqrt(abs(B).^2+K2.^2);
    phi3 = -sqrt(abs(B).^2+K3.^2);
    
    
    cs1=cos(phi1);
    si1=sin(phi1);
    cs2=cos(phi2);
    si2=sin(phi2);
    cs3=cos(phi3);
    si3=sin(phi3);
    
    n1 = real(B)./abs(phi1);
    n2 = imag(B)./abs(phi1);
    n3 = K1./abs(phi1);
    n4 = real(B)./abs(phi2);
    n5 = imag(B)./abs(phi2);
    n6 = K2./abs(phi2);
    n7 = real(B)./abs(phi3);
    n8 = imag(B)./abs(phi3);
    n9 = K3./abs(phi3);
    
    
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
    
    Bd19 = n7.*n7.*(1-cs3)+cs3;
    Bd20 = n7.*n8.*(1-cs3)-n9.*si3;
    Bd21 = n7.*n9.*(1-cs3)+n8.*si3;
    Bd22 = n8.*n7.*(1-cs3)+n9.*si3;
    Bd23 = n8.*n8.*(1-cs3)+cs3;
    Bd24 = n8.*n9.*(1-cs3)-n7.*si3;
    Bd25 = n9.*n7.*(1-cs3)-n8.*si3;
    Bd26 = n9.*n8.*(1-cs3)+n7.*si3;
    Bd27 = n9.*n9.*(1-cs3)+cs3;
    
    
    %Rotation
    %             Mrot=Bd*Mn;
    Mrot(1,:)=Bd1'.*Mt(1,:)+Bd2'.*Mt(2,:)+Bd3'.*Mt(3,:);
    Mrot(2,:)=Bd4'.*Mt(1,:)+Bd5'.*Mt(2,:)+Bd6'.*Mt(3,:);
    Mrot(3,:)=Bd7'.*Mt(1,:)+Bd8'.*Mt(2,:)+Bd9'.*Mt(3,:);
    Mrot(4,:)=Bd10'.*Mt(4,:)+Bd11'.*Mt(5,:)+Bd12'.*Mt(6,:);
    Mrot(5,:)=Bd13'.*Mt(4,:)+Bd14'.*Mt(5,:)+Bd15'.*Mt(6,:);
    Mrot(6,:)=Bd16'.*Mt(4,:)+Bd17'.*Mt(5,:)+Bd18'.*Mt(6,:);
    Mrot(7,:)=Bd19'.*Mt(7,:)+Bd20'.*Mt(8,:)+Bd21'.*Mt(9,:);
    Mrot(8,:)=Bd22'.*Mt(7,:)+Bd23'.*Mt(8,:)+Bd24'.*Mt(9,:);
    Mrot(9,:)=Bd25'.*Mt(7,:)+Bd26'.*Mt(8,:)+Bd27'.*Mt(9,:);
    
    Mrot_n=Mrot+A2b;
    
    C=DC*Mrot_n;
    Mt=D*C;
    
    Mt=Mt-A2b;
    
    
    Mtn=Dx*Mrot;
    
    Mt=Mt+Mtn;
    
    Mn(1,:)=Bd1'.*Mt(1,:)+Bd2'.*Mt(2,:)+Bd3'.*Mt(3,:);
    Mn(2,:)=Bd4'.*Mt(1,:)+Bd5'.*Mt(2,:)+Bd6'.*Mt(3,:);
    Mn(3,:)=Bd7'.*Mt(1,:)+Bd8'.*Mt(2,:)+Bd9'.*Mt(3,:);
    Mn(4,:)=Bd10'.*Mt(4,:)+Bd11'.*Mt(5,:)+Bd12'.*Mt(6,:);
    Mn(5,:)=Bd13'.*Mt(4,:)+Bd14'.*Mt(5,:)+Bd15'.*Mt(6,:);
    Mn(6,:)=Bd16'.*Mt(4,:)+Bd17'.*Mt(5,:)+Bd18'.*Mt(6,:);
    Mn(7,:)=Bd19'.*Mt(7,:)+Bd20'.*Mt(8,:)+Bd21'.*Mt(9,:);
    Mn(8,:)=Bd22'.*Mt(7,:)+Bd23'.*Mt(8,:)+Bd24'.*Mt(9,:);
    Mn(9,:)=Bd25'.*Mt(7,:)+Bd26'.*Mt(8,:)+Bd27'.*Mt(9,:);
    
    Mt=Mn;
    M=Mn;
end

end