P =  peaks(20);
%[x1,y1,z1,c] = sphere3d(-P,pi/2,2*pi,- pi/4,pi/4,8,.4,'surf',.001);
%[x1,y1,z1,c] = sphere3d(-P,pi/2,2*pi,- pi/4,pi/4,8,.4,'mesh',.001);
%[x1,y1,z1,c] = sphere3d(-P,pi/2,2*pi,- pi/4,pi/4,8,.4,'meshc',.001);

[x1,y1,z1,c] = sphere3d(-P,pi/2,2*pi,- pi/4,pi/4,8,.4,'off','spline',.001); 
[x2,y2,z2,c] = sphere3d(-P,pi/2,2*pi,-pi/4,pi/4,15,.4,'off','spline',.001);
[x3,y3,z3,c] = sphere3d(-P,pi/2,2*pi,-pi/4,pi/4,25,.4,'off','spline',.001);
s = surf(x1,y1,z1,c); 
s.EdgeColor = 'none';
hold on
surf(x2,y2,z2,c); hold on
surf(x3,y3,z3,c); 
axis off;
grid off;
set(gca,'DataAspectRatio',[1 1 1])    
set(gca,'XDir','rev','YDir','rev'); 
hold off
