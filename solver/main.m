clear;clc;

N=400;

p=1.0/2;q=1.0/6; r=0;

L=10;
aa=0.02;
CC0=aa^2*1.6e-6;
CC1=2.8e-6;
CC2=1.6e-6; 
EI1=0.0035;
EI2=0.0035;
GIp=0.0035;
EA=1;
MN=2.6e6;
Derta=0.06;
Gama=0.001;
s_int=3*Derta;

[li0,kap1i_0,kap2i_0,taui_0,X0]=shape_free(N,p,q,r,L);
[kap1ii_0,kap2ii_0,tauii_0,XX0]=shape_free2(N,p,q,r,L);
[X,ti,d2i,li]=shape_ini(N,li0,X0);
kap1_i=kap1i_0;
kap2_i=kap2i_0;
tau_i=taui_0;

figure;
plot3(X0(1,1:N+1),X0(2,1:N+1),X0(3,1:N+1),'r-');
hold on
plot3(XX0(1,1:N+1),XX0(2,1:N+1),XX0(3,1:N+1),'b-');
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

Li0=zeros(1,N+1);
for i=1:N
    if i==1||i==N
        Li0(i+1)=Li0(i)+li0(i)/2;
    else
        Li0(i+1)=Li0(i)+li0(i);
    end
end

Tim=0; MaxIter=200; Tol=1e-12; Steps=16000; FinalTim=0.2;
Dtim=FinalTim/Steps; NoSnap=5; SnapInc=Steps/NoSnap;

V=zeros(4*N+3,1);


X_Time_All=zeros(4*N+3,Steps);
Force_ext=zeros(4*N+3,Steps);
Real_V=zeros(4*N+3,Steps);
Time_All=zeros(Steps,1);
PhaseV_All=zeros(Steps,2);
Kapa_All=zeros(N-1,Steps);
Tau_All=zeros(N-1,Steps);
Wm_All=zeros(Steps,6);
Wt_All=zeros(Steps,6);
phi_All=zeros(Steps,N+1);
Ps_All=zeros(Steps,1);
Derta_All=zeros(Steps,1);


for ii=1:Steps

    Tim=Tim+Dtim;
    XP=X;
    MKN=3*MN*Derta;
    MKD=48/(pi^2-4)*MN*Derta;
    H=24*Gama/Derta;
    lam=3*Derta*Gama;

    Ps_All(ii)=s_int;
    Derta_All(ii)=Derta;
    s_int_old=s_int;
    Derta_old=Derta;

    phi=1/2.*(1+tanh((Li0-s_int)./Derta));
    phi_All(i,:)=phi;

    [U,idxf]=BoundaryR(N,Tim,Dtim);

    for jj=1:MaxIter

        [F_extR,K_extR]=Force_extR(N,ti,d2i,li,V,li0,CC1,CC2,CC0,Dtim);

        [F_intR,K_intR]=Force_intR(N,ti,d2i,li,li0,...
            kap1i_0,kap2i_0,taui_0,kap1ii_0,kap2ii_0,tauii_0,...
            p,q,r,phi,EI1,EI2,GIp,EA);
        KR=K_intR+K_extR;
        
        if jj==1
            RHS=-F_extR-F_intR-KR*U;
        else
            RHS=-F_extR-F_intR;
        end

        [LK,UK]=lu(KR(idxf,idxf));
        dU=LK\RHS(idxf);
        dU=UK\dU;

        U(idxf)=dU;
        X=X+U;

        [ti, d2i,li]=Update_everything(N,X,U,ti,d2i,li0);%更新几何
        
        V = (X-XP)/Dtim; 

        if jj==1
            error1=abs(dot(RHS(idxf),dU));
        else
            errormsg=abs(dot(RHS(idxf),dU))/error1;
            if  errormsg<Tol
                break;
            end
        end
        U=U*0;
    end

    [kap1_i,kap2_i,tau_i]=Calculate_curvature(N,p,q,r,ti,d2i,li0);

    for kk=1:MaxIter
        [Fdiss_m,W_T_m,W_M_m,W_Double_m,DFdiss_mm,DFdiss_mt]=Compute_Fdiss(N,kap1_i,kap2_i,tau_i,Li0,s_int,Derta,...
            kap1i_0,kap2i_0,taui_0,kap1ii_0,kap2ii_0,tauii_0,li0,EI1,EI2,GIp,Gama,1);

        [Fdiss_t,W_T_t,W_M_t,W_Double_t,DFdiss_tt,DFdiss_tm]=Compute_Fdiss(N,kap1_i,kap2_i,tau_i,Li0,s_int,Derta,...
            kap1i_0,kap2i_0,taui_0,kap1ii_0,kap2ii_0,tauii_0,li0,EI1,EI2,GIp,Gama,2);

        Rm=s_int-s_int_old-MKN*Fdiss_m*Dtim;
        Rt=Derta-Derta_old-MKD*Fdiss_t*Dtim;

        Kss=1-Dtim*MKN*DFdiss_mm;
        Ksd=-Dtim*MKN/Derta*Fdiss_m-Dtim*MKN*DFdiss_mt;
        Kdd=1-Dtim*MKD/Derta*Fdiss_t-Dtim*MKD*DFdiss_tt;
        Kds=-Dtim*MKD*DFdiss_tm;

        det=Kss*Kdd-Ksd*Kds;
        ds=(-Rm*Kdd+Rt*Ksd)/det;
        dD=(-Rt*Kss+Rm*Kds)/det;

        s_int=s_int+ds;
        Derta=Derta+dD;

        if abs(ds)/Derta+abs(dD)/Derta<Tol&&abs(Rm)/Derta+abs(Rt)/Derta<Tol
            break;
        end
    end

    PhaseV_m=(s_int-s_int_old)/Dtim;
    PhaseV_t=(Derta-Derta_old)/Dtim;
    PhaseV_All(ii,:)=[PhaseV_m PhaseV_t];
    Wm_All(ii,:)=[Fdiss_m,W_T_m,W_M_m,W_Double_m,DFdiss_mm,DFdiss_mt];
    Wt_All(ii,:)=[Fdiss_t,W_T_t,W_M_t,W_Double_t,DFdiss_tt,DFdiss_tm];

    if s_int<3*Derta
        s_int=3*Derta;
    end
    
    X_Time_All(:,ii)=X;
    Force_ext(:,ii)=F_extR;
    Real_V(:,ii)=V;
    Time_All(ii)=Tim;
    Kapa_All(:,ii)=sqrt(kap1_i.^2+kap2_i.^2);
    Tau_All(:,ii)=tau_i;

    if jj==MaxIter
        ii
        jj
        disp(['Iter ', num2str(jj), ': norm(RHS)=', num2str(norm(RHS(idxf))), ...
      ', norm(dU)=', num2str(norm(dU))]);
        error('Newton-Raphson iterations did not converge!');
    end
end





