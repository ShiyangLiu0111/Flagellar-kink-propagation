function [F_int, KK]=Force_intR(N,ti,d2i,li,li_0,...
            kap1i_0,kap2i_0,taui_0,kap1ii_0,kap2ii_0,tauii_0,...
            p,q,r,phi,EI1,EI2,GIp,EA)

F_int=zeros(4*N+3,1);
KK=zeros(4*N+3,4*N+3);

Dui_p0=[zeros(3,3) zeros(3,1) -eye(3) zeros(3,1) eye(3)];
Dui_m0=[-eye(3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)];


Dui_mStem=[-eye(3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)];
Dui_pStem=[zeros(3,3) zeros(3,1) -eye(3) zeros(3,1) eye(3)];
T1=eye(11)+[-Dui_mStem;zeros(8,11)];
TN=eye(11)+[zeros(8,11);Dui_pStem];

phii_p=[zeros(1,3) 0 zeros(1,3) 1 zeros(1,3)];
phii_m=[zeros(1,3) 1 zeros(1,3) 0 zeros(1,3)];

kap1_i=zeros(N-1,1);
kap2_i=zeros(N-1,1);
tau_i =zeros(N-1,1);

Li0=0;

for i=2:N

    li0_p=li_0(i); 
    li0_m=li_0(i-1);
    li0=(li0_p+li0_m)/2;
    
    ti_p=ti(:,i);
    ti_m=ti(:,i-1);
    
    d2i_p=d2i(:,i);
    d2i_m=d2i(:,i-1);
    
    d1i_p=cross(d2i_p,ti_p);
    d1i_m=cross(d2i_m,ti_m);

    La=p+q*(dot(ti_m,ti_p)+dot(d1i_m,d1i_p)+dot(d2i_m,d2i_p));

    kap1_iM=-(1/La+r)*dot(ti_p,d2i_m)/li0;
    kap1_iP=(1/La+r)*dot(ti_m,d2i_p)/li0;

    kap2_iM=(1/La+r)*dot(ti_p,d1i_m)/li0;
    kap2_iP=-(1/La+r)*dot(ti_m,d1i_p)/li0;

    tau_iM=(1/La+r)*dot(d1i_p,d2i_m)/li0;
    tau_iP=-(1/La+r)*dot(d1i_m,d2i_p)/li0;

    kap1_i(i-1)=(kap1_iM+kap1_iP)/2;
    kap2_i(i-1)=(kap2_iM+kap2_iP)/2;
    tau_i(i-1) =(tau_iM+tau_iP)/2;
end

for i=2:N
    dof11_i=dof11(i,1:11);

    li0_p=li_0(i);
    li0_m=li_0(i-1);
    li0=(li0_p+li0_m)/2;

    Li0=Li0+li0;


    Dui_p=Dui_p0/li0_p;
    Dui_m=Dui_m0/li0_m;
    
    li_p=li(i);
    li_m=li(i-1);
    
    ti_p=ti(:,i);
    ti_m=ti(:,i-1);
    
    d2i_p=d2i(:,i);
    d2i_m=d2i(:,i-1);
    
    d1i_p=cross(d2i_p,ti_p);
    d1i_m=cross(d2i_m,ti_m);
    
    kap1_iM=-dot(ti_p,d2i_m)/li0;
    kap1_iP=dot(ti_m,d2i_p)/li0;
    kap1_iSt=(kap1_iP+kap1_iM)/2;

    kap2_iM=dot(ti_p,d1i_m)/li0;
    kap2_iP=-dot(ti_m,d1i_p)/li0;
    kap2_iSt=(kap2_iP+kap2_iM)/2;

    tau_iM=dot(d1i_p,d2i_m)/li0;
    tau_iP=-dot(d1i_m,d2i_p)/li0;
    tau_iSt=(tau_iP+tau_iM)/2;
    
    La=p+q*(dot(ti_m,ti_p)+dot(d1i_m,d1i_p)+dot(d2i_m,d2i_p));
    
    kap1_i0=kap1i_0(i-1);
    kap2_i0=kap2i_0(i-1);
    tau_i0=taui_0(i-1);

    kap1_ii0=kap1ii_0(i-1);
    kap2_ii0=kap2ii_0(i-1);
    tau_ii0=tauii_0(i-1);
    
    Pti_p=eye(3)-ti_p*ti_p';
    Pti_m=eye(3)-ti_m*ti_m';
    
    VarIvLa=-1/La^2*q*((dot(d2i_m,d1i_p)-dot(d1i_m,d2i_p))*phii_m ...
        +(dot(d2i_p,d1i_m)-dot(d1i_p,d2i_m))*phii_p ...
        +(ti_p'*Pti_m-dot(d1i_p,ti_m)*d1i_m'-dot(d2i_p,ti_m)*d2i_m')*Dui_m/li_m ...
        +(ti_m'*Pti_p-dot(d1i_m,ti_p)*d1i_p'-dot(d2i_m,ti_p)*d2i_p')*Dui_p/li_p);
    
    VarVarIvLa=Dui_m'/li_m*(Pti_m*Pti_p+dot(ti_m,ti_p)*(d1i_m*d1i_p'+d2i_m*d2i_p'))*Dui_p/li_p ...
        +(dot(d2i_m,d2i_p)+dot(d1i_m,d1i_p))*(phii_m'*phii_p)...
        +Dui_p'/li_p*(-dot(d1i_m,ti_p)*d2i_p+dot(d2i_m,ti_p)*d1i_p)*(phii_p-phii_m)...
        +Dui_m'/li_m*(-dot(d1i_p,ti_m)*d2i_m+dot(d2i_p,ti_m)*d1i_m)*(phii_m-phii_p);
    VarVarIvLa=VarVarIvLa+VarVarIvLa';
    
    VarVarIvLa=VarVarIvLa ...
        +1/li_m^2*Dui_m'*(-ti_p*ti_m'-ti_m*ti_p'-ti_m'*ti_p*(eye(3)-3*(ti_m*ti_m')) ...
        +dot(d1i_p,ti_m)*(2*d1i_m*ti_m'+ti_m*d1i_m') ...
        +dot(d2i_p,ti_m)*(2*d2i_m*ti_m'+ti_m*d2i_m')...
        -(d1i_m*d1i_p'+d2i_m*d2i_p'))*Dui_m...
        +1/li_p^2*Dui_p'*(-ti_p*ti_m'-ti_m*ti_p'-ti_m'*ti_p*(eye(3)-3*(ti_p*ti_p')) ...
        +dot(d1i_m,ti_p)*(2*d1i_p*ti_p'+ti_p*d1i_p') ...
        +dot(d2i_m,ti_p)*(2*d2i_p*ti_p'+ti_p*d2i_p')...
        -(d1i_p*d1i_m'+d2i_p*d2i_m'))*Dui_p...
        +(-dot(d1i_m,d1i_p)-dot(d2i_m,d2i_p))*(phii_m'*phii_m+phii_p'*phii_p);
    VarVarIvLa=q*VarVarIvLa;
    VarVarIvLa=-1/La^2*VarVarIvLa+2*La*(VarIvLa'*VarIvLa);
    
    %%%%%%%%%%%%%%%%%%%% 1-bending  %%%%%%%%%%%%%%%%%%%%%%%%%%
    B1_p1=-dot(ti_m,ti_p)*d2i_p'*Dui_p/li_p;
    B1_p2=d2i_p'*Pti_m*Dui_m/li_m;
    B1_p3=-d1i_p'*ti_m*phii_p;
    
    B1_m1=-dot(ti_p,ti_m)*d2i_m'*Dui_m/li_m;
    B1_m2=d2i_m'*Pti_p*Dui_p/li_p;
    B1_m3=-d1i_m'*ti_p*phii_m;
    
    B1_M=-(B1_m1+B1_m2+B1_m3)/li0;
    B1_P=(B1_p1+B1_p2+B1_p3)/li0;
    B1_St=(B1_P+B1_M)/2;
    B1M=B1_M*(1/La+r)+kap1_iM*VarIvLa;
    B1P=B1_P*(1/La+r)+kap1_iP*VarIvLa;
    B1=(B1M+B1P)/2;
    
    KG1=Dui_p'/li_p*(dot(ti_p,d2i_m)*(eye(3)-3*(ti_p*ti_p'))+ti_p*d2i_m'+d2i_m*ti_p'...
        +dot(ti_m,ti_p)*(ti_p*d2i_p'+2*d2i_p*ti_p')-d2i_p*ti_m')*Dui_p/li_p;
    
    KG1=KG1-dot(ti_m,d2i_p)*(phii_p'*phii_p);
    KGtemp=-Dui_p'/li_p*(d2i_p*ti_p')*Pti_m*Dui_m/li_m;
    KGtemp=KGtemp+dot(ti_p,ti_m)*Dui_p'/li_p*d1i_p*phii_p;
    KGtemp=KGtemp-Dui_m'/li_m*Pti_m*d1i_p*phii_p;
    KG1=KG1+KGtemp+KGtemp';
    KG1_P=KG1/li0;
    
    KG1=Dui_m'/li_m*(dot(ti_m,d2i_p)*(eye(3)-3*(ti_m*ti_m'))+ti_m*d2i_p'+d2i_p*ti_m'...
        +dot(ti_p,ti_m)*(ti_m*d2i_m'+2*d2i_m*ti_m')-d2i_m*ti_p')*Dui_m/li_m;
    
    KG1=KG1-dot(ti_p,d2i_m)*(phii_m'*phii_m);
    KGtemp=-Dui_m'/li_m*(d2i_m*ti_m')*Pti_p*Dui_p/li_p;
    KGtemp=KGtemp+dot(ti_m,ti_p)*Dui_m'/li_m*d1i_m*phii_m;
    KGtemp=KGtemp-Dui_p'/li_p*Pti_p*d1i_m*phii_m;
    KG1=KG1+KGtemp+KGtemp';
    KG1_M=-KG1/li0;
    KG1=(KG1_M+KG1_P)/2;
    
    KGex=B1_St'*VarIvLa;
    KG1=KG1*(1/La+r)+KGex+KGex'+kap1_iSt*VarVarIvLa;
    
    %%%%%%%%%%%%%%%%%%%% 2-bending  %%%%%%%%%%%%%%%%%%%%%%%%%%
    B2_p1=dot(ti_m,ti_p)*d1i_p'*Dui_p/li_p;
    B2_p2=-d1i_p'*Pti_m*Dui_m/li_m;
    B2_p3=-d2i_p'*ti_m*phii_p;
    
    B2_m1=dot(ti_p,ti_m)*d1i_m'*Dui_m/li_m;
    B2_m2=-d1i_m'*Pti_p*Dui_p/li_p;
    B2_m3=-d2i_m'*ti_p*phii_m;
    
    B2_M=-(B2_m1+B2_m2+B2_m3)/li0;
    B2_P=(B2_p1+B2_p2+B2_p3)/li0;
    B2_St=(B2_P+B2_M)/2;
    
    B2M=B2_M*(1/La+r)+kap2_iM*VarIvLa;
    B2P=B2_P*(1/La+r)+kap2_iP*VarIvLa;
    B2=B2_St*(1/La+r)+kap2_iSt*VarIvLa;
    
    KG2=Dui_p'/li_p*(-dot(ti_p,d1i_m)*(eye(3)-3*(ti_p*ti_p'))-ti_p*d1i_m'-d1i_m*ti_p'...
        +dot(ti_m,ti_p)*(-ti_p*d1i_p'-2*d1i_p*ti_p')+d1i_p*ti_m')*Dui_p/li_p;
    
    KG2=KG2+dot(ti_m,d1i_p)*(phii_p'*phii_p);
    KGtemp=Dui_p'/li_p*(d1i_p*ti_p')*Pti_m*Dui_m/li_m;
    KGtemp=KGtemp+dot(ti_p,ti_m)*Dui_p'/li_p*d2i_p*phii_p;
    KGtemp=KGtemp-Dui_m'/li_m*Pti_m*d2i_p*phii_p;
    KG2=KG2+KGtemp+KGtemp';
    KG2_P=KG2/li0;
    
    KG2=Dui_m'/li_m*(-dot(ti_m,d1i_p)*(eye(3)-3*(ti_m*ti_m'))-ti_m*d1i_p'-d1i_p*ti_m'...
        +dot(ti_p,ti_m)*(-ti_m*d1i_m'-2*d1i_m*ti_m')+d1i_m*ti_p')*Dui_m/li_m;
    
    KG2=KG2+dot(ti_p,d1i_m)*(phii_m'*phii_m);
    KGtemp=Dui_m'/li_m*(d1i_m*ti_m')*Pti_p*Dui_p/li_p;
    KGtemp=KGtemp+dot(ti_m,ti_p)*Dui_m'/li_m*d2i_m*phii_m;
    KGtemp=KGtemp-Dui_p'/li_p*Pti_p*d2i_m*phii_m;
    KG2=KG2+KGtemp+KGtemp';
    KG2_M=-KG2/li0;
    KG2=(KG2_M+KG2_P)/2;
    
    KGex=B2_St'*VarIvLa;
    KG2=KG2*(1/La+r)+KGex+KGex'+kap2_iSt*VarVarIvLa;
    
    % %%%%%%%%%%%%%%%%%%%% 3-twisting  %%%%%%%%%%%%%%%%%%%%%%%%%%
    B3_p1=dot(d1i_m,d1i_p)*phii_p;
    B3_p2=-dot(d2i_m,d2i_p)*phii_m;
    B3_p3=dot(d2i_p,ti_m)*d1i_m'*Dui_m/li_m;
    B3_p4=dot(d1i_m,ti_p)*d2i_p'*Dui_p/li_p;
    
    B3_m1=dot(d1i_p,d1i_m)*phii_m;
    B3_m2=-dot(d2i_p,d2i_m)*phii_p;
    B3_m3=dot(d2i_m,ti_p)*d1i_p'*Dui_p/li_p;
    B3_m4=dot(d1i_p,ti_m)*d2i_m'*Dui_m/li_m;
    
    B3_M=-(B3_m1+B3_m2+B3_m3+B3_m4)/li0;
    B3_P=(B3_p1+B3_p2+B3_p3+B3_p4)/li0;
    B3_St=(B3_P+B3_M)/2;

    B3M=B3_M*(1/La+r)+tau_iM*VarIvLa;
    B3P=B3_P*(1/La+r)+tau_iP*VarIvLa;
    B3=B3_St*(1/La+r)+tau_iSt*VarIvLa;
    
    KG3=Dui_p'*(-dot(ti_p,d1i_m)*(ti_p*d2i_p'+2*d2i_p*ti_p')+d2i_p*d1i_m'...
        +dot(ti_p,d2i_m)*(ti_p*d1i_p'+2*d1i_p*ti_p')-d1i_p*d2i_m')*Dui_p/li_p^2;
    
    KG3=KG3+dot(d1i_m,d2i_p)*(phii_p'*phii_p+phii_m'*phii_m);
    KGtemp=-ti_p'*ti_m/li_p/li_m*(Dui_p'*d2i_p*d1i_m'*Dui_m);
    KGtemp=KGtemp-Dui_p'/li_p*(dot(ti_p,d2i_m)*d2i_p+dot(ti_p,d1i_m)*d1i_p)*(phii_p-phii_m);
    KGtemp=KGtemp+dot(d1i_p,d2i_m)*(phii_p'*phii_m);
    KG3=KG3+KGtemp+KGtemp';
    KG3_P=KG3/li0;
    
    KG3=Dui_m'*(-dot(ti_m,d1i_p)*(ti_m*d2i_m'+2*d2i_m*ti_m')+d2i_m*d1i_p'...
        +dot(ti_m,d2i_p)*(ti_m*d1i_m'+2*d1i_m*ti_m')-d1i_m*d2i_p')*Dui_m/li_m^2;
    
    KG3=KG3+dot(d1i_p,d2i_m)*(phii_m'*phii_m+phii_p'*phii_p);
    KGtemp=-ti_m'*ti_p/li_m/li_p*(Dui_m'*d2i_m*d1i_p'*Dui_p);
    KGtemp=KGtemp-Dui_m'/li_m*(dot(ti_m,d2i_p)*d2i_m+dot(ti_m,d1i_p)*d1i_m)*(phii_m-phii_p);
    KGtemp=KGtemp+dot(d1i_m,d2i_p)*(phii_m'*phii_p);
    KG3=KG3+KGtemp+KGtemp';
    KG3_M=-KG3/li0;
    KG3=(KG3_M+KG3_P)/2;
    
    KGex=B3_St'*VarIvLa;
    KG3=KG3*(1/La+r)+KGex+KGex'+tau_iSt*VarVarIvLa;
    
    % %%%%%%%%%%%%%%%%%%%% streching  %%%%%%%%%%%%%%%%%%%%%%%%%%
    Bs_p=li_p*ti_p'*Dui_p;
    Ni_p=EA*(li_p^2-1)/2;
    Bs_m=li_m*ti_m'*Dui_m;
    Ni_m=EA*(li_m^2-1)/2;
    
    F_ints=(Bs_p'*Ni_p*li0_p+Bs_m'*Ni_m*li0_m)/2;
    
    KMs=EA/2*(Bs_p'*Bs_p*li0_p+Bs_m'*Bs_m*li0_m);
    
    KGs=((Dui_p'*Dui_p)*Ni_p*li0_p+(Dui_m'*Dui_m)*Ni_m*li0_m)/2;
    
    %%%%%%%%%%%%%%%%%%%%    Phase Transition   %%%%%%%%%%%%%%%%%%%%%%%%%%
    Mb1i=EI1*(kap1_iSt*(1/La+r)-kap1_ii0-phi(i)*(kap1_i0-kap1_ii0));
    F_intb1=(B1'*Mb1i)*li0;
    KM1=EI1*li0*(B1'*B1);
    KG1=KG1*Mb1i*li0;

    Mb2i=EI2*(kap2_iSt*(1/La+r)-kap2_ii0-phi(i)*(kap2_i0-kap2_ii0));
    F_intb2=(B2'*Mb2i)*li0;
    KM2=EI2*li0*(B2'*B2);
    KG2=KG2*Mb2i*li0;

    Mti=GIp*(tau_iSt*(1/La+r)-tau_ii0-phi(i)*(tau_i0-tau_ii0));
    F_intt=(B3'*Mti)*li0;
    KM3=GIp*li0*(B3'*B3);
    KG3=KG3*Mti*li0;
    %%%%%%%%%%%%%%%%%%%%     sum   %%%%%%%%%%%%%%%%%%%%%%%%%%

    K=KM1+KM2+KM3+KMs+KG1+KG2+KG3+KGs;
    F_int_=F_intb1+F_intb2+F_intt+F_ints;
    
    if i==2
        KK(dof11_i,dof11_i)=KK(dof11_i,dof11_i)+T1'*K*T1;
        F_int(dof11_i)=F_int(dof11_i)+T1'*F_int_;
     elseif i==N
        KK(dof11_i,dof11_i)=KK(dof11_i,dof11_i)+TN'*K*TN;
        F_int(dof11_i)=F_int(dof11_i)+TN'*F_int_;
    else
        KK(dof11_i,dof11_i)=KK(dof11_i,dof11_i)+K;
        F_int(dof11_i)=F_int(dof11_i)+F_int_;
    end  
end

end