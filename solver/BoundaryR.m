function [URp,idxf]=BoundaryR(N,Tim,Dtim)
twistv=100;
if Tim<=0.01
    BC=[1  0
      2  0
      3  0
      4  2*pi*twistv*Dtim
      5  0
      6  0
      7  0];
else
  BC=[1  0
      2  0
      3  0
      4  -2*pi*twistv*Dtim
      5  0
      6  0
      7  0];
end
  [URp,idxf]=Form_up(N,BC);
end
function [URp,idxf]=Form_up(N,bc)
URp=zeros(4*N+3,1);
URp(bc(:,1))=bc(:,2);
idxf=1:4*N+3;
idxf(bc(:,1))=[];
end
