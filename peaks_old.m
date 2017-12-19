M1=dlmread('solids/overlap1.dat',' '); x=M1(:,1);y=M1(:,2);z=M1(:,3);
I1=(abs(y(:))<1e-12); I2=(abs(x(:))<1e-12);
[max(z(I1)), max(z(I2))]

C = textread('chosen.solution.txt', '%s','delimiter', ' '); 
xlab=['[',C{10},C{11},C{12},']']; ylab=['[',C{13},C{14},C{15},']'];
tri=delaunay(x,y); 
f1=figure('Position', [100, 100, 640, 480]);
trisurf(tri,x,y,z,'LineStyle','none');
ax=gca; ax.FontSize=15;
ax.XTick = [-0.5,-0.25,0,0.25,0.5]; ax.YTick = [-0.5,-0.25,0,0.25,0.5];
ax.XLabel.String = xlab; ax.YLabel.String = ylab; ax.ZLabel.String = 'overlap volume (A^3)'; 
ax.XLabel.Rotation = 15; ax.YLabel.Rotation = -20;
ax.ZLabel.Position(1)=ax.ZLabel.Position(1)-0.02;
%ax.Title.String = 'Plane 1'; 
print(f1,'OverlapPlot0','-dpng','-r100');
close(f1);
