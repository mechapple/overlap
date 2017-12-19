M1=dlmread('solids/overlap.dat',' '); x=M1(:,1);y=M1(:,2);z=M1(:,3);
X=[x;x+1;x;x+1];Y=[y;y;y+1;y+1];Z=[z;z;z;z]; R2=unique([X Y Z],'rows');
TOL=1e-6;
R1=R2;
%R1=R2(R2(:,1)>=-TOL & R2(:,1)<=(1+TOL) & R2(:,2)>=-TOL & R2(:,2)<=(1+TOL),:);
C = textread('chosen.solution.txt', '%s','delimiter', ' ');
xlab=['x=[',C{10},C{11},C{12},']']; ylab=['y=[',C{13},C{14},C{15},']'];
tri=delaunay(R1(:,1),R1(:,2));

f1=figure('Position', [100, 100, 640, 480]);
%trisurf(tri,R1(:,1),R1(:,2),R1(:,3),'LineStyle','none');
ax=gca; ax.FontSize=15;
ax.XTick = min(R1(:,1)):0.25:max(R1(:,1)); ax.YTick = min(R1(:,2)):0.25:max(R1(:,2));
ax.XLabel.String = xlab; ax.YLabel.String = ylab; ax.ZLabel.String = 'overlap volume (Å^3)';
%ax.XLabel.Rotation = 15; ax.YLabel.Rotation = -20;
%ax.ZLabel.Position(1)=ax.ZLabel.Position(1)-0.02;
%ax.Title.String = 'Plane 1';
x=R1(:,1);y=R1(:,2);z=R1(:,3);
ti=-0.5:0.01:1.5;
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(x,y,z,XI,YI,'cubic');
tricontour([x y],tri,z,20)
print(f1,'OverlapPlot','-dpng','-r100');
close(f1);

f2=figure('Position', [100, 100, 640, 480], 'DefaultAxesFontSize',18);
t1=-0.5:0.0001:1.5;
I=t1(:)>=0.0 & t1(:)<=1.0;
w1=interp2(XI,YI,ZI,t1,0); [val1, idx1] = max(w1'.*I);
w2=interp2(XI,YI,ZI,0,t1); [val2, idx2] = max(w2.*I);
w3=interp2(XI,YI,ZI,t1,t1); [val3, idx3] = max(w3'.*I);
maxlocs=[t1(idx1) val1 val1/t1(idx1);t1(idx2) val2 val2/t1(idx2);t1(idx3) val3 val3/t1(idx3) ];
dlmwrite('max_overlap.dat',maxlocs,'delimiter',' ','precision','%.6f');
plot(t1,w1,'r',t1,w2,'g',t1,w3,'b');
legend('y=0','x=0','x=y');
xlabel('t'); ylabel('overlap volume (Å^3)');
print(f2,'Slices','-dpng','-r100');
close(f2);
