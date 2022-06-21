close all
clc
clear all; 


xmin = 0; 
xmax = 10; 
ymin = 0; 
ymax = 10; 
ds = 0.5; 
order = 1; 
offset = xmax-0.2*ds; %xmax-ds/2-1.25; %xmax-ds*2; 

Lx = xmax-xmin; 
Ly = ymax-ymin; 
nx = Lx/ds; 
ny = Ly/ds; 

nelem = nx*ny*2; 
nnodes = (nx+1)*(ny+1);
coords = zeros(nnodes,2); 

fid1 = fopen('grid1.dat','w'); 
fprintf(fid1,'%i %i %i\n',nnodes,nelem,order); 

fid2 = fopen('grid2.dat','w'); 
fprintf(fid2,'%i %i %i\n',nnodes,nelem,order); 

% form nodes
n = 0; 
x = xmin:ds:xmax; 
y = ymin:ds:ymax; 
for i = 1:(ny+1)
    for j = 1:(nx+1)
        n=n+1; 
        coords(n,1) = x(j); 
        coords(n,2) = y(i); 
       
        fprintf(fid1,'%.8e %.8e\n',x(j),y(i));
        fprintf(fid2,'%.8e %.8e\n',x(j)+offset,y(i));
    end
end

np1 = n; % number of p=1 nodes


% form elements
for i = 1:ny % row
    for j = 1:nx % col
        a = (i-1)*(nx+1)+j;
        b = a+1; 
        c = i*(nx+1)+j;
        d = c+1; 
        
        fprintf(fid1,'%i %i %i\n',a,b,c); 
        fprintf(fid1,'%i %i %i\n',d,c,b); 
        fprintf(fid2,'%i %i %i\n',a,d,c); 
        %fprintf(fid2,'%i %i %i\n',a,b,c); 
        %fprintf(fid2,'%i %i %i\n',d,c,b); 
        fprintf(fid2,'%i %i %i\n',a,b,d); 

        xtmp1 = [coords(a,1) coords(b,1) coords(c,1) coords(a,1)]; 
        ytmp1 = [coords(a,2) coords(b,2) coords(c,2) coords(a,2)];         
        xtmp2 = [coords(d,1) coords(c,1) coords(b,1) coords(d,1)]; 
        ytmp2 = [coords(d,2) coords(c,2) coords(b,2) coords(d,2)];  

%         figure(1)
%         plot(xtmp1,ytmp1,'k',xtmp2,ytmp2,'k');
%         hold on; 
    end
end

% boundary nodes
nbnodes=0;
for i= 1:nnodes
  if(coords(i,1)==xmin || coords(i,1) == xmax || coords(i,2)==ymin || coords(i,2) == ymax)     nbnodes=nbnodes+1;
  end
end
fprintf(fid1,"%i\n",nbnodes);
fprintf(fid2,"%i\n",nbnodes);
for i = 1:nnodes
    if(coords(i,1)==xmin || coords(i,1) == xmax || coords(i,2)==ymin || coords(i,2) == ymax)
        fprintf(fid1,'%i 1\n',i); 
        fprintf(fid2,'%i 1\n',i); 
    end
end
fclose(fid1); 
fclose(fid2); 
