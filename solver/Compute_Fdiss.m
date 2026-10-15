function [Fdiss,W_T,W_M,W_Double,DFdiss1,DFdiss2]=Compute_Fdiss(N,kap1_i,kap2_i,tau_i,Li0,s_int,Derta,...
        kap1i_0,kap2i_0,taui_0,kap1ii_0,kap2ii_0,tauii_0,li_0,EI1,EI2,GIp,Gama,XX)

Li_0=0;
W_T=0;
W_M=0;
W_Double=0;
DFdiss_M1=0;
DFdiss_T1=0;
DW_Double1=0;
DFdiss_M2=0;
DFdiss_T2=0;
DW_Double2=0;

gp=[-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, ...
    0.405845151377397, 0.741531185599394, 0.949107912342759];
gw=[0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, ...
    0.381830050505119, 0.279705391489277, 0.129484966168870];
n_points = 7;

for i=2:N

    li0_p=li_0(i); 
    li0_m=li_0(i-1);
    li0=(li0_p+li0_m)/2;
    J=li0/2;


    s = zeros(1,n_points);
    sint = zeros(1,n_points);
    phi = zeros(1,n_points);
    Dphi = zeros(1,n_points);
    DDphi = zeros(1,n_points);
    DDDphi = zeros(1,n_points);
    
    for j=1:n_points
        s(j)=Li_0+(gp(j)+1)/2*li0;
        sint(j)=(s(j)-s_int)/Derta;
        phi(j)=0.5*(1+tanh(sint(j)));
        Dphi(j)=(2/Derta)*phi(j)*(1-phi(j));
        DDphi(j)=(4/Derta^2)*phi(j)*(1-phi(j))*(1-2*phi(j));
        DDDphi(j)=(8/Derta^3)*phi(j)*(1-phi(j))*(1-6*phi(j)+6*phi(j)^2);
    end


    kap1 = zeros(1,n_points);
    kap2 = zeros(1,n_points);
    tau = zeros(1,n_points);
    kap1_0 = zeros(1,n_points);
    kap1_00 = zeros(1,n_points);
    kap2_0 = zeros(1,n_points);
    kap2_00 = zeros(1,n_points);
    tau_0 = zeros(1,n_points);
    tau_00 = zeros(1,n_points);

    if i-1==1
        for j=1:n_points
            if s(j) < Li0(i)
                w1=(Li0(i)-s(j))/li0_p;
                kap1(j)=kap1_i(i-1)-w1*(kap1_i(i)-kap1_i(i-1));
                kap2(j)=kap2_i(i-1)-w1*(kap2_i(i)-kap2_i(i-1));
                tau(j)=tau_i(i-1)-w1*(tau_i(i)-tau_i(i-1));
                kap1_0(j)=kap1i_0(i-1)-w1*(kap1i_0(i)-kap1i_0(i-1));
                kap1_00(j)=kap1ii_0(i-1)-w1*(kap1ii_0(i)-kap1ii_0(i-1));
                kap2_0(j)=kap2i_0(i-1)-w1*(kap2i_0(i)-kap2i_0(i-1));
                kap2_00(j)=kap2ii_0(i-1)-w1*(kap2ii_0(i)-kap2ii_0(i-1));
                tau_0(j)=taui_0(i-1)-w1*(taui_0(i)-taui_0(i-1));
                tau_00(j)=tauii_0(i-1)-w1*(tauii_0(i)-tauii_0(i-1));
            elseif s(j) == Li0(i)
                kap1(j)=kap1_i(i-1);
                kap2(j)=kap2_i(i-1);
                tau(j)=tau_i(i-1);
                kap1_0(j)=kap1i_0(i-1);
                kap1_00(j)=kap1ii_0(i-1);
                kap2_0(j)=kap2i_0(i-1);
                kap2_00(j)=kap2ii_0(i-1);
                tau_0(j)=taui_0(i-1);
                tau_00(j)=tauii_0(i-1);
            else 
                w3=(s(j)-Li0(i))/li0_p;
                kap1(j)=w3*kap1_i(i)+(1-w3)*kap1_i(i-1);
                kap2(j)=w3*kap2_i(i)+(1-w3)*kap2_i(i-1);
                tau(j)=w3*tau_i(i)+(1-w3)*tau_i(i-1);
                kap1_0(j)=w3*kap1i_0(i)+(1-w3)*kap1i_0(i-1);
                kap1_00(j)=w3*kap1ii_0(i)+(1-w3)*kap1ii_0(i-1);
                kap2_0(j)=w3*kap2i_0(i)+(1-w3)*kap2i_0(i-1);
                kap2_00(j)=w3*kap2ii_0(i)+(1-w3)*kap2ii_0(i-1);
                tau_0(j)=w3*taui_0(i)+(1-w3)*taui_0(i-1);
                tau_00(j)=w3*tauii_0(i)+(1-w3)*tauii_0(i-1);
            end
        end
        
    elseif i-1==N-1
        for j=1:n_points
            if s(j) < Li0(i)
                w1=(Li0(i)-s(j))/li0_m;
                kap1(j)=w1*kap1_i(i-2)+(1-w1)*kap1_i(i-1);
                kap2(j)=w1*kap2_i(i-2)+(1-w1)*kap2_i(i-1);
                tau(j)=w1*tau_i(i-2)+(1-w1)*tau_i(i-1);
                kap1_0(j)=w1*kap1i_0(i-2)+(1-w1)*kap1i_0(i-1);
                kap1_00(j)=w1*kap1ii_0(i-2)+(1-w1)*kap1ii_0(i-1);
                kap2_0(j)=w1*kap2i_0(i-2)+(1-w1)*kap2i_0(i-1);
                kap2_00(j)=w1*kap2ii_0(i-2)+(1-w1)*kap2ii_0(i-1);
                tau_0(j)=w1*taui_0(i-2)+(1-w1)*taui_0(i-1);
                tau_00(j)=w1*tauii_0(i-2)+(1-w1)*tauii_0(i-1);
            elseif s(j) == Li0(i) 
                kap1(j)=kap1_i(i-1);
                kap2(j)=kap2_i(i-1);
                tau(j)=tau_i(i-1);
                kap1_0(j)=kap1i_0(i-1);
                kap1_00(j)=kap1ii_0(i-1);
                kap2_0(j)=kap2i_0(i-1);
                kap2_00(j)=kap2ii_0(i-1);
                tau_0(j)=taui_0(i-1);
                tau_00(j)=tauii_0(i-1);
            else
                w3=(s(j)-Li0(i))/li0_m;
                kap1(j)=kap1_i(i-1)+w3*(kap1_i(i-1)-kap1_i(i-2));
                kap2(j)=kap2_i(i-1)+w3*(kap2_i(i-1)-kap2_i(i-2));
                tau(j)=tau_i(i-1)+w3*(tau_i(i-1)-tau_i(i-2));
                kap1_0(j)=kap1i_0(i-1)+w3*(kap1i_0(i-1)-kap1i_0(i-2));
                kap1_00(j)=kap1ii_0(i-1)+w3*(kap1ii_0(i-1)-kap1ii_0(i-2));
                kap2_0(j)=kap2i_0(i-1)+w3*(kap2i_0(i-1)-kap2i_0(i-2));
                kap2_00(j)=kap2ii_0(i-1)+w3*(kap2ii_0(i-1)-kap2ii_0(i-2));
                tau_0(j)=taui_0(i-1)+w3*(taui_0(i-1)-taui_0(i-2));
                tau_00(j)=tauii_0(i-1)+w3*(tauii_0(i-1)-tauii_0(i-2));
            end
        end
        
    else
        for j=1:n_points
            if s(j) < Li0(i)
                w1=(Li0(i)-s(j))/li0_m;
                kap1(j)=w1*kap1_i(i-2)+(1-w1)*kap1_i(i-1);
                kap2(j)=w1*kap2_i(i-2)+(1-w1)*kap2_i(i-1);
                tau(j)=w1*tau_i(i-2)+(1-w1)*tau_i(i-1);
                kap1_0(j)=w1*kap1i_0(i-2)+(1-w1)*kap1i_0(i-1);
                kap1_00(j)=w1*kap1ii_0(i-2)+(1-w1)*kap1ii_0(i-1);
                kap2_0(j)=w1*kap2i_0(i-2)+(1-w1)*kap2i_0(i-1);
                kap2_00(j)=w1*kap2ii_0(i-2)+(1-w1)*kap2ii_0(i-1);
                tau_0(j)=w1*taui_0(i-2)+(1-w1)*taui_0(i-1);
                tau_00(j)=w1*tauii_0(i-2)+(1-w1)*tauii_0(i-1);
            elseif s(j) == Li0(i)
                kap1(j)=kap1_i(i-1);
                kap2(j)=kap2_i(i-1);
                tau(j)=tau_i(i-1);
                kap1_0(j)=kap1i_0(i-1);
                kap1_00(j)=kap1ii_0(i-1);
                kap2_0(j)=kap2i_0(i-1);
                kap2_00(j)=kap2ii_0(i-1);
                tau_0(j)=taui_0(i-1);
                tau_00(j)=tauii_0(i-1);
            else
                w3=(s(j)-Li0(i))/li0_p;
                kap1(j)=w3*kap1_i(i)+(1-w3)*kap1_i(i-1);
                kap2(j)=w3*kap2_i(i)+(1-w3)*kap2_i(i-1);
                tau(j)=w3*tau_i(i)+(1-w3)*tau_i(i-1);
                kap1_0(j)=w3*kap1i_0(i)+(1-w3)*kap1i_0(i-1);
                kap1_00(j)=w3*kap1ii_0(i)+(1-w3)*kap1ii_0(i-1);
                kap2_0(j)=w3*kap2i_0(i)+(1-w3)*kap2i_0(i-1);
                kap2_00(j)=w3*kap2ii_0(i)+(1-w3)*kap2ii_0(i-1);
                tau_0(j)=w3*taui_0(i)+(1-w3)*taui_0(i-1);
                tau_00(j)=w3*tauii_0(i)+(1-w3)*tauii_0(i-1);
            end
        end
    end

    kap101 = zeros(1,n_points);
    kap102 = zeros(1,n_points);
    kap103 = zeros(1,n_points);
    kap201 = zeros(1,n_points);
    kap202 = zeros(1,n_points);
    kap203 = zeros(1,n_points);
    tau01 = zeros(1,n_points);
    tau02 = zeros(1,n_points);
    tau03 = zeros(1,n_points);
    
    for j=1:n_points
        kap101(j)=kap1_00(j)+phi(j)*(kap1_0(j)-kap1_00(j));
        kap102(j)=kap1_00(j)+phi(j)*(kap1_0(j)-kap1_00(j));
        kap103(j)=kap1_00(j)+phi(j)*(kap1_0(j)-kap1_00(j));
        kap201(j)=kap2_00(j)+phi(j)*(kap2_0(j)-kap2_00(j));
        kap202(j)=kap2_00(j)+phi(j)*(kap2_0(j)-kap2_00(j));
        kap203(j)=kap2_00(j)+phi(j)*(kap2_0(j)-kap2_00(j));
        tau01(j)=tau_00(j)+phi(j)*(tau_0(j)-tau_00(j));
        tau02(j)=tau_00(j)+phi(j)*(tau_0(j)-tau_00(j));
        tau03(j)=tau_00(j)+phi(j)*(tau_0(j)-tau_00(j));
    end

    Li_0=Li_0+li0;

    if XX==1
        for j=1:n_points
            W_T=W_T-GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*Dphi(j)*gw(j)*J;
            W_M=W_M-EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*Dphi(j)*gw(j)*J-...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*Dphi(j)*gw(j)*J;
            W_Double=W_Double+6*Gama*Derta*Dphi(j)*DDphi(j)*gw(j)*J;

            DFdiss_M1=DFdiss_M1-EI1*(kap1_0(j)-kap1_00(j))^2*Dphi(j)^2*gw(j)*J+...
                EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*DDphi(j)*gw(j)*J-...
                EI2*(kap2_0(j)-kap2_00(j))^2*Dphi(j)^2*gw(j)*J+...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*DDphi(j)*gw(j)*J;

            DFdiss_T1=DFdiss_T1-GIp*(tau_0(j)-tau_00(j))^2*Dphi(j)^2*gw(j)*J+...
                GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*DDphi(j)*gw(j)*J;

            DW_Double1=DW_Double1-6*Gama*Derta*DDphi(j)^2*gw(j)*J-...
                6*Gama*Derta*Dphi(j)*DDDphi(j)*gw(j)*J;

            DFdiss_M2=DFdiss_M2-EI1*(kap1_0(j)-kap1_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J-...
                EI2*(kap2_0(j)-kap2_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J;

            DFdiss_T2=DFdiss_T2-GIp*(tau_0(j)-tau_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J;

            DW_Double2=DW_Double2-12*Gama*Dphi(j)*DDphi(j)*gw(j)*J-...
                6*Gama*Derta*DDphi(j)^2*sint(j)*gw(j)*J-...
                6*Gama*Derta*Dphi(j)*DDDphi(j)*sint(j)*gw(j)*J;
        end
        
    elseif XX==2
        for j=1:n_points
            W_T=W_T-GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*Dphi(j)*gw(j)*J*sint(j);
            W_M=W_M-EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*Dphi(j)*gw(j)*J*sint(j)-...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*Dphi(j)*gw(j)*J*sint(j);
            W_Double=W_Double+3*Gama*Dphi(j)^2*gw(j)*J+...
                6*Gama*Derta*Dphi(j)*DDphi(j)*sint(j)*gw(j)*J;

            DFdiss_M1=DFdiss_M1-EI1*(kap1_0(j)-kap1_00(j))^2*Dphi(j)^2*sint(j)^2*gw(j)*J+...
                EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*(2*sint(j)/Derta*Dphi(j)+sint(j)^2*DDphi(j))*gw(j)*J-...
                EI2*(kap2_0(j)-kap2_00(j))^2*Dphi(j)^2*sint(j)^2*gw(j)*J+...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*(2*sint(j)/Derta*Dphi(j)+sint(j)^2*DDphi(j))*gw(j)*J;

            DFdiss_T1=DFdiss_T1-GIp*(tau_0(j)-tau_00(j))^2*Dphi(j)^2*sint(j)^2*gw(j)*J+...
                GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*(2*sint(j)/Derta*Dphi(j)+sint(j)^2*DDphi(j))*gw(j)*J;

            DW_Double1=DW_Double1-6*Gama/Derta*Dphi(j)^2*gw(j)*J-...
                24*Gama*Dphi(j)*DDphi(j)*sint(j)*gw(j)*J-...
                6*Gama*Derta*DDphi(j)^2*sint(j)^2*gw(j)*J-...
                6*Gama*Derta*Dphi(j)*DDDphi(j)*sint(j)^2*gw(j)*J;

            DFdiss_M2=DFdiss_M2-EI1*(kap1_0(j)-kap1_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                EI1*(kap1_0(j)-kap1_00(j))*(kap1(j)-kap101(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J-...
                EI2*(kap2_0(j)-kap2_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                EI2*(kap2_0(j)-kap2_00(j))*(kap2(j)-kap201(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J;

            DFdiss_T2=DFdiss_T2-GIp*(tau_0(j)-tau_00(j))^2*Dphi(j)^2*sint(j)*gw(j)*J+...
                GIp*(tau_0(j)-tau_00(j))*(tau(j)-tau01(j))*(Dphi(j)/Derta+sint(j)*DDphi(j))*gw(j)*J;

            DW_Double2=DW_Double2-12*Gama*Dphi(j)*DDphi(j)*gw(j)*J-...
                6*Gama*Derta*DDphi(j)^2*sint(j)*gw(j)*J-...
                6*Gama*Derta*Dphi(j)*DDDphi(j)*sint(j)*gw(j)*J;
        end
    end

end

DFdiss1=DFdiss_M1+DFdiss_T1+DW_Double1;
DFdiss2=DFdiss_M2+DFdiss_T2+DW_Double2;
Fdiss=W_T+W_M+W_Double;

end