%{ 
mu*du/dt=u*( t - u);
u(-1)=3; -1<=t<=2
%}
clear; clc; close all;
%mu = 0.01; M=280 - ломается
%mu=0.1; M=87; 
%mu=0.05; M=905; %904
mu=0.05; M=118; %117-118; 50
a_=-1;
b=2;
tau=(b-a_)/(M - 1); %временной шаг

U1=zeros(M,1); 
U2=zeros(M,1);
T=zeros(M,1);

T(1)=-1;
for m=2:M
   T(m)=T(1)+m*tau;
end

method = 1; % 1 - CROS1; 2 - CROS2

if (method == 1) %CROS1
    a=(1+1i)/2;
    Un=3;
    for m=1:M
        U1(m)=Un;
        w=(-1/mu*U1(m)*(U1(m)-(T(m)+tau/2))) / (1+tau*a*(1/mu*(2*U1(m)-T(m))));
        Un=Un+tau*real(w);
    end
end

method = 2;

if ( method == 2) %CROS2
    a11 = 0.1 + (sqrt(11)/30)*1i;
    a22 = 0.2 + 0.1i;
    b1 = 0.1941430241155180 - 0.2246898944678803i;
    b2 = 0.8058569758844820 - 0.8870089521907592i;
    c21 = 0.2554708972958462 - 0.2026195833570109i;
    a21 = 0.5617645150714754 - 1.148223341045841i;
        
    w1=zeros(2,1);
    w2=zeros(2,1);
    Um(1,:)=[3 -1];
    
    syms u1 u2;
    u=[u1 u2];
    F=[-1/mu*u1*(u1-u2); 1];
    Fu=jacobian(F,u);
    
    f1=sym(Fu);
    f2=sym(F);
    
    v1=symvar(Fu);
    v2=symvar(F);
    
    g1=@(X1) double(subs(f1,v1,X1));
    g2=@(X2) double(subs(f2,v2,X2));
   
    for m=2:M
        G1=g1([Um(m-1,1) Um(m-1,2)]);
        G2=g2([Um(m-1,1) Um(m-1,2)]);
        w1=( eye(2) - a11*tau*G1 )\G2;
        
        e=Um(m-1,1)+tau*real(a21*w1(1)');
        b=Um(m-1,2)+tau*real(a21*w1(2)');
        c=Um(m-1,1)+tau*real(c21*w1(1)');
        d=Um(m-1,2)+tau*real(c21*w1(2)');
        
        G1s=g1([e b]);
        G2s=g2([c d]);
        w2=( eye(2) - a22*tau*G1s ) \G2s ;
        
        Um(m,:)=Um(m-1,:)+tau*real( b1*w1'+b2*w2') ;
        
    end
    U2=Um(:,1);
    T=Um(:,2);
end

%аналитическое решение задачи
K = 1000;
U = zeros(K,1);

T2=zeros(K,1);
T2(1)=-1;
tau2=(b-a_)/(K - 1);
for m=2:K
   T2(m)=T2(1)+m*tau2;
end  
for m=1:K
    if T2(m) == -1
        U(m,1) = 3;
    end
    if ((T2(m) > -1) && (T2(m) <=1))
        U(m,1) = 0;
    end
    if T2(m) > 1
        U(m,1) = T2(m);
    end
end

figure(method);
cros1 = plot(T,U1,'g');
hold on;
cros2 = plot(T,U2);
anal = plot(T2,U,'r');
set(cros1,'LineWidth',2);
set(cros2,'LineWidth',2);
set(anal,'LineWidth',2);
legend('CROS1','CROS2','Аналитическое решение')

