function [dofnum]=dof_i(i,j)
[II,JJ]=meshgrid((i-1)*4,j);
dofnum=II+JJ;
dofnum=dofnum(:)';
end 