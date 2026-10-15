filename='X_Time_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'zone T="Frame"\r\n');
    fprintf(fid,'STRANDID=%d\r\n',1);
    fprintf(fid,'I=%d\r\n',N+1);
    fprintf(fid,'ZONETYPE=Ordered\r\n');
    fprintf(fid,'DATAPACKING=POINT\r\n');
    fprintf(fid,'SOLUTIONTIME=%.5f\r\n',(i-1)/Steps*Tim);
    for j=1:N+1
        fprintf(fid,'%d %d %d %d\n',X_Time_All(dof_i(j,1),i),X_Time_All(dof_i(j,2),i),X_Time_All(dof_i(j,3),i),phi_All(i,j));
    end
end
fclose(fid);

filename='Force_ext.txt';
fid=fopen(filename,'w');
for i=1:Steps
    for j=1:4*N+3
        fprintf(fid,'%d\n',Force_ext(j,i));
    end
end
fclose(fid);

filename='Real_V.txt';
fid=fopen(filename,'w');
for i=1:Steps
    for j=1:4*N+3
        fprintf(fid,'%d\n',Real_V(j,i));
    end
end
fclose(fid);

filename='Time_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d\n',Time_All(i,1));
end
fclose(fid);

filename='s_int_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d\n',Ps_All(i,1));
end
fclose(fid);

filename='Derta_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d\n',Derta_All(i,1));
end
fclose(fid);

filename='PhaseV_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d %d\n',PhaseV_All(i,1),PhaseV_All(i,2));
end
fclose(fid);

filename='Kapa_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    for j=1:N-1
        fprintf(fid,'%d\n',Kapa_All(j,i));
    end
end
fclose(fid);

filename='Tau_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    for j=1:N-1
        fprintf(fid,'%d\n',Tau_All(j,i));
    end
end
fclose(fid);


filename='Wm_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d %d %d %d %d %d\n',Wm_All(i,:));
end
fclose(fid);

filename='Wt_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'%d %d %d %d %d %d\n',Wt_All(i,:));
end
fclose(fid);

%%
for i=1:Steps
    for j=1:N+1
        X_Time_All(j+407*(i-1),4)=phi(i,j);
    end
end
%%
filename='X_Time_All.txt';
fid=fopen(filename,'w');
for i=1:Steps
    fprintf(fid,'zone T="Frame"\r\n');
    fprintf(fid,'STRANDID=%d\r\n',1);
    fprintf(fid,'I=%d\r\n',N+1);
    fprintf(fid,'ZONETYPE=Ordered\r\n');
    fprintf(fid,'DATAPACKING=POINT\r\n');
    fprintf(fid,'SOLUTIONTIME=%.5f\r\n',(i-1)/Steps*Tim);
    for j=1:N+1
        fprintf(fid,'%d %d %d %d\n',X_Time_All(j+407*(i-1),1),X_Time_All(j+407*(i-1),2),X_Time_All(j+407*(i-1),3),X_Time_All(j+407*(i-1),4));
    end
end
fclose(fid);