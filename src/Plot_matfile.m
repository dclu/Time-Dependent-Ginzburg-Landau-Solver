ddt=100
dx=0.2;
dy=0.2;

ina=1;
ind=18;
% filename=strcat('MYGLaaaatt',num2str(ina));
% % filename=strcat('MYGLaaaatt',num2str(ina));
% filename=strcat(filename,num2str(ind*ddt));
filename='F:\Eralop\Phy Projects\GL Vortex\matlab file\MyGL\lattice\Article Pic\vortex_lattice_final';
% filename=strcat('F:\Eralop\Phy Projects\GL Vortex\matlab file\MyGL\lattice\mat\MYGLapp',num2str(ina));
% filename=strcat(filename,num2str(ind));
s=load(filename);


NNx = 2;
NNy = 4;

[Nxp,Nyp]=size(s.x);
Nx=Nxp-2;
Ny=Nyp-2;

psip = zeros(Nx*NNx,Ny*NNy);
psip2 = zeros(Nx*NNx,Ny*NNy);
psip3 = zeros(Nx*NNx,Ny*NNy);
psipt = zeros(Nx*NNx,Ny*NNy);
pB = zeros(Nx*NNx,Ny*NNy);
pjx = zeros(Nx*NNx,Ny*NNy);
pjy = zeros(Nx*NNx,Ny*NNy);





for ix =1:NNx
for jy = 1:NNy
    for idx = 1:Nx
        for jdy = 1:Ny
            psip(Nx*(ix-1)+idx,Ny*(jy-1)+jdy)=abs(s.x(idx+1,jdy+1))^2;
            psip2(Nx*(ix-1)+idx,Ny*(jy-1)+jdy)=abs(s.y(idx+1,jdy+1))^2;
            psip3(Nx*(ix-1)+idx,Ny*(jy-1)+jdy)=abs(s.z(idx+1,jdy+1))^2;
            
            pB(Nx*(ix-1)+idx,Ny*(jy-1)+jdy) = real( (s.n(idx+2,jdy+1)-s.n(idx+1-1,jdy+1))/2/dx - (s.l(idx+1,jdy+2)-s.l(idx+1,jdy+1-1))/2/dy);
            %pjx(Nx*(ix-1)+idx,Ny*(jy-1)+jdy) = (ax(idx+1,jdy+2)+ax(idx+1,jdy+0)-2*ax(idx+1,jdy+1))/dx^2- (ay(idx+1+1,jdy+1+1)+ay(idx+1-1,jdy+1-1)-ay(idx+1-1,jdy+1+1)-ay(idx+1+1,jdy+1-1))/dx/dy/4 ;
            %pjy(Nx*(ix-1)+idx,Ny*(jy-1)+jdy) = (ay(idx+1+1,jdy+1)+ay(idx+1-1,jdy+1)-2*ay(idx+1,jdy+1))/dy^2- (ax(idx+1+1,jdy+1+1)+ax(idx+1-1,jdy+1-1)-ax(idx+1-1,jdy+1+1)-ax(idx+1+1,jdy+1-1))/dx/dy/4 ;
        end
    end 
end
end
% fre=0;
% for idx = 2:Nx+1
%     for jdy = 2:Ny+1
%         fre = fre-((a1^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1)+(a1^2*ax(idx,jdy)^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1+(a1^2*ay(idx,jdy)^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1+(sqrt(a1*a2*b1*b2)*g*psi2(idx,jdy)*conj(psi(idx,jdy)))/(b1*b2)+(a1^2*psi(idx,jdy)^2*conj(psi(idx,jdy))^2)/(2*b1)+(sqrt(a1*a2*b1*b2)*g*psi(idx,jdy)*conj(psi2(idx,jdy)))/(b1*b2)-(a2^2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/b2+(a1*a2*m1*ax(idx,jdy)^2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(b2*m2)+(a1*a2*m1*ay(idx,jdy)^2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(b2*m2)+(a2^2*psi2(idx,jdy)^2*conj(psi2(idx,jdy))^2)/(2*b2)+(a1^2*(ax(idx,jdy+1)-ax(idx,jdy))/dy^2)/b1+(1i*a1^2*ay(idx,jdy)*conj(psi(idx,jdy))*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(b1*k1)+(1i*a1*a2*m1*ay(idx,jdy)*conj(psi2(idx,jdy))*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(b2*k1*m2)-(1i*a1^2*ay(idx,jdy)*psi(idx,jdy)*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1)+(a1^2*(psi(idx,jdy+1)-psi(idx,jdy))/dy*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1^2)-(1i*a1*a2*m1*ay(idx,jdy)*psi2(idx,jdy)*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b2*k1*m2)+(a1*a2*m1*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b2*k1^2*m2)-(2*a1^2*(ax(idx,jdy+1)-ax(idx,jdy))/dy*(ay(idx+1,jdy)-ay(idx,jdy))/dx)/b1+(a1^2*(ay(idx+1,jdy)-ay(idx,jdy))/dx^2)/b1+(1i*a1^2*ax(idx,jdy)*conj(psi(idx,jdy))*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(b1*k1)+(1i*a1*a2*m1*ax(idx,jdy)*conj(psi2(idx,jdy))*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(b2*k1*m2)-(1i*a1^2*ax(idx,jdy)*psi(idx,jdy)*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1)+(a1^2*(psi(idx+1,jdy)-psi(idx,jdy))/dx*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1^2)-(1i*a1*a2*m1*ax(idx,jdy)*psi2(idx,jdy)*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b2*k1*m2)+(a1*a2*m1*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b2*k1^2*m2);
%     end
% end
% fre






LL=4;


rat=1.67;

Ly = LL/rat;
Lx = LL*rat;

ppt=transpose(psip3);

x_vecp = linspace(0,NNx*Lx,NNx*Nx);
y_vecp = linspace(0,NNy*Ly,NNy*Ny);
[xx,yy] = meshgrid(x_vecp,y_vecp);

%%%begin细化
x_f = linspace(0,NNx*Lx,NNx*Nx*8);
y_f = linspace(0,NNy*Ly,NNy*Ny*8);
[xfine,yfine] = meshgrid(x_f,y_f);

ppf=interp2(xx,yy,ppt,xfine,yfine,'cubic');
%%%end细化


% nx=4;
% ny=4;
% 
% x_c = linspace(0,NNx*Lx*8-2*nx*dx,NNx*Nx*8-2*nx);
% y_c = linspace(0,NNy*Ly*8-2*ny*dy,NNy*Ny*8-2*ny);
% [xc,yc] = meshgrid(x_c,y_c);
% 
% pc = zeros(Nx*NNx*8-2*nx,Ny*NNy*8-2*ny);
% for i = 1:Nx*NNx*8-2*nx
%     for j = 1:Ny*NNy*8-2*ny
%         pc(i,j)=ppf(i+nx,j+ny);
%         
%     end
% 
% end 


[c,h]=contourf(xfine,yfine,ppf,4);

set(gca,'DataAspectRatio',[1 1 1])
colormap(parula(7));
xlabel('L_x/\lambda');
ylabel('L_y/\lambda');


% mesh(xfine,yfine,ppf);
% view(-15,17)
% zlabel('|\psi_s|^2');
% %zlabel('\phi^2');
% 
% xlabel('L_x/\lambda');
% ylabel('L_y/\lambda');
% set(gca,'DataAspectRatio',[1 1 0.10])
% colormap(jet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%meshc
% handles = meshc(xfine,yfine,ppf);
% % handles = meshc(xc,yc,pc);
% 
% %view(0,0)
% view(-15,17)
% 
% zlabel('|\psi_s|^2');
% %zlabel('\phi^2');
% 
% xlabel('L_x/\lambda');
% ylabel('L_y/\lambda');
% set(gca,'DataAspectRatio',[1 1 0.1])
% 
% 
% alpha(1.0);
% %%%颜色
% colormap(jet);
% 
% %%%end颜色
%  
% % handles is a 2-element array of handles: the surface plot and the contours
% hContour = handles(2); % get the handle to the contour lines
% hContour.ContourZLevel = 0.2; % set the contour's Z position (default: hAxes.ZLim(1)=-10)
%  
% % We can also customize other aspects of the contour lines, for example:
% hContour.LineWidth = 0.4; % set the contour lines' width (default: 0.5)


