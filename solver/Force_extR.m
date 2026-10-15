function [F_ext,K_ext]=Force_extR(N,ti,d2i,li,V,li_0,CC1,CC2,CC0,Dtim)
F_ext=zeros(4*N+3,1);

K_ext=zeros(4*N+3,4*N+3);

Dui_p0=[zeros(3,3) zeros(3,1) -eye(3) zeros(3,1) eye(3)];
Dui_m0=[-eye(3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)];


Fui_m=[eye(3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)]/2;
Fui=[zeros(3,3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)];
Fui_p=[zeros(3,3) zeros(3,1) eye(3) zeros(3,1) eye(3)]/2;

Fphii_m=[zeros(1,3) 1 zeros(1,3) 0 zeros(1,3)];
Fphii=[zeros(1,3) 1 zeros(1,3) 1 zeros(1,3)]/2;
Fphii_p=[zeros(1,3) 0 zeros(1,3) 1 zeros(1,3)];

Eui_m1=(1/Dtim)*[eye(3) zeros(3,8)];
Eui=(1/Dtim)*[zeros(3,4) eye(3) zeros(3,4)];
Eui_p1=(1/Dtim)*[zeros(3,8) eye(3)];

Dui_mStem=[-eye(3) zeros(3,1) eye(3) zeros(3,1) zeros(3,3)];
Dui_pStem=[zeros(3,3) zeros(3,1) -eye(3) zeros(3,1) eye(3)];
T1=eye(11)+[-Dui_mStem;zeros(8,11)];
TN=eye(11)+[zeros(8,11);Dui_pStem];

for i=2:N

    dof11_i=dof11(i,1:11);
    Vi=V(dof11_i,1);

    li_m=li(i-1);
    li_p=li(i);

    li0_p=li_0(i);
    li0_m=li_0(i-1);
    li0=(li0_p+li0_m)/2;

    Dui_p=Dui_p0/li0_p;
    Dui_m=Dui_m0/li0_m;

    ti_m=ti(:,i-1);
    ti_p=ti(:,i);
    ti_mp=(ti_m+ti_p)/2;

    d2i_m=d2i(:,i-1);
    d2i_p=d2i(:,i);
    d2i_mp=(d2i_p+d2i_m)/2;
    
    d1i_p=cross(d2i_p,ti_p);
    d1i_m=cross(d2i_m,ti_m);
    d1i_mp=(d1i_p+d1i_m)/2;

    Pti_m=eye(3)-ti_m*ti_m'; 
    Pti_p=eye(3)-ti_p*ti_p';

    
    Vi_m=(Vi(1:3,1)+Vi(5:7,1))/2;
    Vi_p=(Vi(5:7,1)+Vi(9:11,1))/2;
    Vi_mp=Vi(5:7,1);

    Fi_m=CC1*dot(Vi_m,d1i_m)*d1i_m...
        +CC1*dot(Vi_m,d2i_m)*d2i_m...
        +CC2*dot(Vi_m,ti_m)*ti_m;
    Fi_p=CC1*dot(Vi_p,d1i_p)*d1i_p...
        +CC1*dot(Vi_p,d2i_p)*d2i_p...
        +CC2*dot(Vi_p,ti_p)*ti_p;
    Fi=CC1*dot(Vi_mp,d1i_mp)*d1i_mp...
        +CC1*dot(Vi_mp,d2i_mp)*d2i_mp...
        +CC2*dot(Vi_mp,ti_mp)*ti_mp;
    Mi_m=CC0*Vi(4);
    Mi_p=CC0*Vi(8);
    Mi=(Mi_m+Mi_p)/2;
    Fe=li0/12*(Fui_m'*Fi_m+4*Fui'*Fi+Fui_p'*Fi_p...
        +Mi_m*Fphii_m'+4*Mi*Fphii'+Mi_p*Fphii_p');%11*1
    

    %刚度
    Vartvti_m=ti_m*( Vi_m'*Pti_m*(Dui_m/li_m) + ti_m'*(0.5*(Eui_m1+Eui)) ) ...
             +dot(ti_m,Vi_m)*Pti_m*(Dui_m/li_m);%3*11
    F_Var_m=CC1*Vartvti_m+CC2*(0.5*(Eui_m1+Eui)-Vartvti_m); 

    Vartvti_i=ti_mp*( Vi_mp'*(0.5*(Pti_m*(Dui_m/li_m)+Pti_p*(Dui_p/li_p))) + ti_mp'*Eui ) ...
             +dot(ti_mp,Vi_mp)*(0.5*(Pti_m*(Dui_m/li_m)+Pti_p*(Dui_p/li_p)));
    F_Var=CC1*Vartvti_i+CC2*(Eui-Vartvti_i);

    Vartvti_p=ti_p*( Vi_p'*Pti_p*(Dui_p/li_p) + ti_p'*(0.5*(Eui_p1+Eui)) ) ...
             +dot(ti_p,Vi_p)*Pti_p*(Dui_p/li_p);
    F_Var_p=CC1*Vartvti_p+CC2*(0.5*(Eui_p1+Eui)-Vartvti_p); 

    
    F_Var_m=F_Var_m'*Fui_m;
    F_Var_p=F_Var_p'*Fui_p;
    F_Var=F_Var'*Fui;

    Ke=li0/6*(F_Var_m+4*F_Var+F_Var_p...
      +Fphii_m'*CC0*Fphii_m+4*Fphii'*CC0*Fphii+Fphii_p'*CC0*Fphii_p);%11*11

    if i==2
        K_ext(dof11_i,dof11_i)=K_ext(dof11_i,dof11_i)+T1'*Ke*T1;
        F_ext(dof11_i)=F_ext(dof11_i)+T1'*Fe;
     elseif i==N
        K_ext(dof11_i,dof11_i)=K_ext(dof11_i,dof11_i)+TN'*Ke*TN;
        F_ext(dof11_i)=F_ext(dof11_i)+TN'*Fe;
    else
        K_ext(dof11_i,dof11_i)=K_ext(dof11_i,dof11_i)+Ke;
        F_ext(dof11_i)=F_ext(dof11_i)+Fe;
    end  
end

end 
