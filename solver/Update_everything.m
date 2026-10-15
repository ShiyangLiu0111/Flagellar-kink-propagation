function [ti,d2i,li]=Update_everything(N,X,U,tip,d2ip,li0)
ti=zeros(3,N);
Ti=ti;
d2i=ti;
li=zeros(1,N);

for i=1:N

    if i==1||i==N
        Ti(:,i)=2*(X(dof_i(i+1,1:3))-X(dof_i(i,1:3)))/li0(i);
    else
        Ti(:,i)=(X(dof_i(i+1,1:3))-X(dof_i(i,1:3)))/li0(i); 
    end
    
    li(i)=sqrt(dot(Ti(:,i),Ti(:,i)));
    ti(:,i)=Ti(:,i)/li(i);

    d2i(:,i)=para_base(tip(:,i),ti(:,i),d2ip(:,i));
    
    thetai=U(dof_i(i,4));
    d2i(:,i)=rotation(ti(:,i),thetai,d2i(:,i));
    d2i(:,i)=d2i(:,i)/norm(d2i(:,i));
end

end
