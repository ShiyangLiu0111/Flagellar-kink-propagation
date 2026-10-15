function [X,tip,d2ip,li]=shape_ini(N,li0,X0)

  X=X0;
  
  lli=zeros(1,N);
  li=lli;
  tip=zeros(3,N);
  d2ip=zeros(3,N);

  for ii=1:N
    dX=X(:,ii+1)-X(:,ii);
    lli(ii)=sqrt(dot(dX,dX));
    tip(:,ii)=dX/lli(ii);
    li(ii)=lli(ii)/li0(ii);
  end
  li(1)=li(1)*2;
  li(end)=li(end)*2;
  theta=0*X(1,:);
  X=[X;theta];
  X=X(:);
  X(end)=[];
  
  for ii=2:N
      ti_p=tip(:,ii);
      ti_m=tip(:,ii-1);
      
      d2ip_temp=cross(ti_m,ti_p);
      d2ip(:,ii-1)=d2ip_temp/norm(d2ip_temp);
  end
  d2ip(:,N)=para_base(tip(:,N-1),tip(:,N),d2ip(:,N-1));
  d2ip(:,N)=d2ip(:,N)/norm(d2ip(:,N));
  
end
