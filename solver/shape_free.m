function [li0,kap1i_0,kap2i_0,taui_0,X0]=shape_free(N,p,q,r,L)
  K0=2; 

  R0=0.2161;                           
  Pitch=2.1611;                       
  theta=atan(2*pi*R0/Pitch);            

  L0=L.*cos(theta);               
  
  n=N-1;                                      
  seg=(1:2:2*n-1)/(2*n);          
  seg=[0 seg (2*n)/(2*n)];
  seg=seg.*L0;    

  xx= R0*(1-exp(-K0*seg.^2)).*cos((2*pi/Pitch)*seg);
  yy=-R0*(1-exp(-K0*seg.^2)).*sin((2*pi/Pitch)*seg);                
  zz=seg;

  X0=[xx;yy;zz];
  
  ti0=zeros(3,N);
  d2i0=zeros(3,N);

  kap1i_0=-0*ones(N-1,1);  
  kap2i_0=-0*ones(N-1,1);
  taui_0=-0*ones(N-1,1);
  li0=zeros(N,1);
  for ii=1:N
    dX=X0(:,ii+1)-X0(:,ii);
    li0(ii)=sqrt(dot(dX,dX));    
    ti0(:,ii)=dX/li0(ii);    
  end

  li0(1)=2*li0(1);           
  li0(N)=2*li0(N);
  
  for ii=2:N
      ti_p=ti0(:,ii);
      ti_m=ti0(:,ii-1);
      
      d2i0_temp=cross(ti_m,ti_p); 
      d2i0(:,ii-1)=d2i0_temp/norm(d2i0_temp);
  end
  d2i0(:,N)=para_base(ti0(:,N-1),ti0(:,N),d2i0(:,N-1));
  d2i0(:,N)=d2i0(:,N)/norm(d2i0(:,N));
  
  
  for ii=2:N  
    ti_p=ti0(:,ii);
    ti_m=ti0(:,ii-1);
    
    d2i_p=d2i0(:,ii);
    d2i_m=d2i0(:,ii-1);

    d1i_p=cross(d2i_p,ti_p);
    d1i_m=cross(d2i_m,ti_m);

    La=p+q*(dot(ti_m,ti_p)+dot(d1i_m,d1i_p)+dot(d2i_m,d2i_p));

    li0_p=li0(ii); 
    li0_m=li0(ii-1);
    l_i0=(li0_p+li0_m)/2;
    
    kap1i_0(ii-1)=(dot(ti_m,d2i_p)-dot(ti_p,d2i_m))/2/l_i0*(1/La+r);
    
    kap2i_0(ii-1)=(-dot(ti_m,d1i_p)+dot(ti_p,d1i_m))/2/l_i0*(1/La+r);
    
    taui_0(ii-1)=(-dot(d1i_m,d2i_p)+dot(d1i_p,d2i_m))/2/l_i0*(1/La+r);    
  end



end 
