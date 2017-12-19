M=dlmread('Crackface2_points.dat',' ');x=M(:,1);y=M(:,2);
k=boundary(x,y); 
A=polyarea(x(k),y(k))
