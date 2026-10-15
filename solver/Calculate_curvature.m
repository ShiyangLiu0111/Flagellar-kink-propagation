function [kap1_i,kap2_i,tau_i]=Calculate_curvature(N,p,q,r,ti,d2i,li_0)

kap1_i=zeros(N-1,1);
kap2_i=zeros(N-1,1);
tau_i =zeros(N-1,1);

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
end