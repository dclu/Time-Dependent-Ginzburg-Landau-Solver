%%% Corresponding to the file "Formal Procedure_LDC Ver_Nematic order.nb"


clear
clc
close all

% psi(t,x,y) ax(t,x,y) ay(t,x,y)

% d^2 psi/dx^2 = (psi(t,x+dx,y)+psi(t,x-dx,y)-2*psi(t,x,y))/dx^2
% d^2 psi/dy^2 = (psi(t,x,y+dy)+psi(t,x,y-dy)-2*psi(t,x,y))/dy^2
% d^2 ax/dx^2 = (ax(t,x+dx,y)+ax(t,x-dx,y)-2*ax(t,x,y))/dx^2
% d^2 ax/dy^2 = (ax(t,x,y+dy)+ax(t,x,y-dy)-2*ax(t,x,y))/dy^2
% d^2 ay/dx^2 = (ay(t,x+dx,y)+ay(t,x-dx,y)-2*ay(t,x,y))/dx^2
% d^2 ay/dy^2 = (ay(t,x,y+dy)+ay(t,x,y-dy)-2*ay(t,x,y))/dy^2
% d^2 ax/dxdy = (ax(t,x+dx,y+dy)+ax(t,x-dx,y-dy)-ax(t,x-dx,y+dy)-ax(t,x+dx,y-dy))/dx/dy/4
% d^2 ay/dxdy = (ay(t,x+dx,y+dy)+ay(t,x-dx,y-dy)-ay(t,x-dx,y+dy)-ay(t,x+dx,y-dy))/dx/dy/4

% define constants
phi = 2*2*3.1415926;
k1=10;
ss=10;

a1 = 1;
b1 = 1;
m1 = 1;


a2 = 1.;
b2 = 1;
m2 = 2;

g1 = 0; %|s|^2 |d|^2%
g2 = 0;  % sc^2 d^2
g3 = 0.2;% anisotropy dynamic


gg=2;

a3 = 1/gg^2;
b3 = 1/gg^2;
m3 = 1*gg^2;


l1 = 0.;%nematic-sd
l2 = 0;

l3 = 0.3/gg;%nematic-s
l4 = 0.3/gg;%nematic-d
%%% psi complex order psi2: real order
%%%

% Size of unit cell
LL=4;


% Displacement of X and Y mesh
dx = 0.2;
dy = 0.2;


% time from 1-100;
dt = 0.005;
t = 20.;
Nt = floor(t/dt);

% define the order parameter and magnetic field (t,x,y)


% % initial state
% for idx = 2:Nx+1
%     for jdy = 2:Ny+1
%         psi(idx,jdy)=((idx-1)*dx-Lx/2-0.5*dx)/sqrt(((idx-1)*dx-Lx/2-0.5*dx)^2+((jdy-1)*dy-Ly/2-0.5*dx)^2)+1i*((jdy-1)*dy-Ly/2-0.5*dy)/sqrt(((idx-1)*dx-Lx/2-0.5*dx)^2+((jdy-1)*dy-Ly/2-0.5*dy)^2);
%         psi2(idx,jdy)=((idx-1)*dx-Lx/2-0.5*dx)/sqrt(((idx-1)*dx-Lx/2-0.5*dx)^2+((jdy-1)*dy-Ly/2-0.5*dx)^2)+1i*((jdy-1)*dy-Ly/2-0.5*dy)/sqrt(((idx-1)*dx-Lx/2-0.5*dx)^2+((jdy-1)*dy-Ly/2-0.5*dy)^2);
%     end
% end 
%psi(1,floor(idx/2),floor(jdy/2))=0;

% set periodic boundary condition
% for jdy = 2:Ny+1
%     psi(1,Nx+2,jdy) = psi(1,2,jdy)*exp(+1i*jdy*dy*phi/2/Lx);
%     psi(1,1,jdy) = psi(1,Nx+1,jdy)*exp(-1i*jdy*dy*phi/2/Lx);
%     ax(1,Nx+2,jdy) = ax(1,2,jdy);
%     ax(1,1,jdy) = ax(1,Nx+1,jdy);
%     ay(1,Nx+2,jdy) = ay(1,2,jdy) + phi/2/Ly;
%     ay(1,1,jdy) = ay(1,Nx+1,jdy) - phi/2/Ly;
% end
% 
% for idx = 2:Nx+1
%     psi(1,idx,Ny+2) = psi(1,idx,2)*exp(-1i*idx*dx*phi/2/Ly);
%     psi(1,idx,1) = psi(1,idx,Ny+1)*exp(+1i*idx*dx*phi/2/Ly);
%     ax(1,idx,Ny+2) = ax(1,idx,2) - phi/2/Lx;
%     ax(1,idx,1) = ax(1,idx,Ny+1) + phi/2/Lx;
%     ay(1,idx,Ny+2) = ay(1,idx,2);
%     ay(1,idx,1) = ay(1,idx,Ny+1);
% end

%%%change the lattice
Na = 1;
% spa=linspace(1.0,2.1,Na);
%%% take snapshot begin
%interval
ddt=100;
pp = [];
pp2 = [];
pp3 = [];
pax = [];
pay = [];
%%% take snapshot end


%change the lattice

%define the energy
EEE2=[];
EEE3=[];
% tt0=clock;
for lda = 1: Na

dt = 0.005;
rat = 1.67;
Ly = LL / rat;
Lx = LL * rat;

Nx=floor(Lx/dx);
Ny=floor(Ly/dy);
    
    
%one step
psi = zeros(Nx+2, Ny+2);
psi2 = zeros(Nx+2, Ny+2);
psi3 = zeros(Nx+2, Ny+2);
ax = zeros(Nx+2, Ny+2);
ay = zeros(Nx+2, Ny+2);

%next step
psik = zeros(Nx+2, Ny+2);
psi2k = zeros(Nx+2, Ny+2);
psi3k = zeros(Nx+2, Ny+2);
axk = zeros(Nx+2, Ny+2);
ayk = zeros(Nx+2, Ny+2);    
    

dpsi = 0;
dpsi2 = 0;
dpsi3 = 0;
dax = 0;
day = 0;
    
    

% init1=-sqrt(2)/2+sqrt(2)*randn(Nx,Ny);
% init2=-sqrt(2)/2+sqrt(2)*randn(Nx,Ny);
% init3=-sqrt(2)/2+sqrt(2)*randn(Nx,Ny);
% init4=-sqrt(2)/2+sqrt(2)*randn(Nx,Ny);
% init5=-1+2*randn(Nx,Ny);

% initial state
for idx = 2:Nx+1
    for jdy = 2:Ny+1
        init=((idx-1)*dx-Lx/2-0.1*dx)/sqrt(((idx-1)*dx-Lx/2-0.1*dx)^2+((jdy-1)*dy-Ly/2-0.1*dx)^2)+1i*((jdy-1)*dy-Ly/2-0.1*dy)/sqrt(((idx-1)*dx-Lx/2-0.1*dx)^2+((jdy-1)*dy-Ly/2-0.1*dy)^2);
        psi(idx,jdy)=0.2*init;
        psi2(idx,jdy)=init;
        psi3(idx,jdy)=0.5; %
    end
end





%%% integrate this

err=1;
errk=10;
count=0;
while(err>1e-9)
    
%    if mod(kdt,ddt) == 0
% %        pp=[pp, psi];
% %        pp2= [pp2, psi2];
% %        pp3= [pp3, psi3];
% %        pax= [pax,ax];
% %        pay= [pay,ay];
%        Nt*Na-(kdt+(lda-1)*Nt)
%        [psi(10,10),psi2(10,10),psi3(10,10)]


   
   

    %%set periodic boundary
    for jdy = 2:Ny+1
        psi(Nx+2,jdy) = psi(2,jdy)*exp(+1i*(jdy-1)*dy*phi/2/Ly);
        psi(1,jdy) = psi(Nx+1,jdy)*exp(-1i*(jdy-1)*dy*phi/2/Ly);
        psi2(Nx+2,jdy) = psi2(2,jdy)*exp(+1i*(jdy-1)*dy*phi/2/Ly);
        psi2(1,jdy) = psi2(Nx+1,jdy)*exp(-1i*(jdy-1)*dy*phi/2/Ly);
        psi3(Nx+2,jdy) = psi3(2,jdy);
        psi3(1,jdy) = psi3(Nx+1,jdy);
        
        ax(Nx+2,jdy) = ax(2,jdy);
        ax(1,jdy) = ax(Nx+1,jdy);
        ay(Nx+2,jdy) = ay(2,jdy) + phi/2/Ly/k1;
        ay(1,jdy) = ay(Nx+1,jdy) - phi/2/Ly/k1;
    end

    for idx = 2:Nx+1
        psi(idx,Ny+2) = psi(idx,2)*exp(-1i*(idx-1)*dx*phi/2/Lx);
        psi(idx,1) = psi(idx,Ny+1)*exp(+1i*(idx-1)*dx*phi/2/Lx);
        psi2(idx,Ny+2) = psi2(idx,2)*exp(-1i*(idx-1)*dx*phi/2/Lx);
        psi2(idx,1) = psi2(idx,Ny+1)*exp(+1i*(idx-1)*dx*phi/2/Lx);
        psi3(idx,Ny+2) = psi3(idx,2);
        psi3(idx,1) = psi3(idx,Ny+1);
        
        ax(idx,Ny+2) = ax(idx,2) - phi/2/Lx/k1;
        ax(idx,1) = ax(idx,Ny+1) + phi/2/Lx/k1;
        ay(idx,Ny+2) = ay(idx,2);
        ay(idx,1) = ay(idx,Ny+1);
    end
    ax(Nx+2,Ny+2) = ax(2,2)-phi/2/Lx/k1;
    ax(1,1) = ax(Nx+1,Ny+1)+phi/2/Lx/k1;
    ay(Nx+2,Ny+2) = ay(2,2)+phi/2/Ly/k1;
    ay(1,1) = ay(Nx+1,Ny+1)-phi/2/Ly/k1;
    
    fre = 0;
    err=0;
    %integrate it
    for idx = 2:Nx+1
        for jdy = 2:Ny+1
            
            dpsi = (psi(idx,jdy)-ax(idx,jdy)^2*psi(idx,jdy)-ay(idx,jdy)^2*psi(idx,jdy)-2*g3*m1*ax(idx,jdy)^2*psi2(idx,jdy)+2*g3*m1*ay(idx,jdy)^2*psi2(idx,jdy)-(sqrt(a3/b3)*l1*psi2(idx,jdy)*psi3(idx,jdy))/a1-(1i*sqrt(a3/b3)*l2*psi2(idx,jdy)*psi3(idx,jdy))/a1-(a3*l3*psi(idx,jdy)*psi3(idx,jdy)^2)/(a1*b3)-psi(idx,jdy)^2*conj(psi(idx,jdy))-(2*g2*psi2(idx,jdy)^2*conj(psi(idx,jdy)))/b1-(g1*psi(idx,jdy)*psi2(idx,jdy)*conj(psi2(idx,jdy)))/b1-(1i*psi(idx,jdy)*(ay(idx,jdy+1)-ay(idx,jdy))/dy)/k1+(2*1i*g3*m1*psi2(idx,jdy)*(ay(idx,jdy+1)-ay(idx,jdy))/dy)/k1-(2*1i*ay(idx,jdy)*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/k1+(4*1i*g3*m1*ay(idx,jdy)*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/k1+(psi(idx,jdy+1)+psi(idx,jdy-1)-2*psi(idx,jdy))/dy^2/k1^2-(2*g3*m1*(psi2(idx,jdy+1)+psi2(idx,jdy-1)-2*psi2(idx,jdy))/dy^2)/k1^2-(1i*psi(idx,jdy)*(ax(idx+1,jdy)-ax(idx,jdy))/dx)/k1-(2*1i*g3*m1*psi2(idx,jdy)*(ax(idx+1,jdy)-ax(idx,jdy))/dx)/k1-(2*1i*ax(idx,jdy)*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/k1-(4*1i*g3*m1*ax(idx,jdy)*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/k1+(psi(idx+1,jdy)+psi(idx-1,jdy)-2*psi(idx,jdy))/dx^2/k1^2+(2*g3*m1*(psi2(idx+1,jdy)+psi2(idx-1,jdy)-2*psi2(idx,jdy))/dx^2)/k1^2)*dt;
            psik(idx,jdy) = psi(idx,jdy) + dpsi;
            
            dpsi2 = (-((2*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*ax(idx,jdy)^2*psi(idx,jdy))/(a2^2*b1))+(2*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*ay(idx,jdy)^2*psi(idx,jdy))/(a2^2*b1)+(sqrt(a2*b1)*sqrt(a1*b2)*psi2(idx,jdy))/(a2*b1)-(a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*ax(idx,jdy)^2*psi2(idx,jdy))/(a2^2*b1*m2)-(a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*ay(idx,jdy)^2*psi2(idx,jdy))/(a2^2*b1*m2)-(sqrt(a2*b1)*sqrt(a1*b2)*sqrt(a3/b3)*l1*psi(idx,jdy)*psi3(idx,jdy))/(a2^2*b1)+(1i*sqrt(a2*b1)*sqrt(a1*b2)*sqrt(a3/b3)*l2*psi(idx,jdy)*psi3(idx,jdy))/(a2^2*b1)-(a3*sqrt(a2*b1)*sqrt(a1*b2)*l4*psi2(idx,jdy)*psi3(idx,jdy)^2)/(a2^2*b1*b3)-(a1*sqrt(a2*b1)*sqrt(a1*b2)*g1*psi(idx,jdy)*psi2(idx,jdy)*conj(psi(idx,jdy)))/(a2^2*b1^2)-(2*a1*sqrt(a2*b1)*sqrt(a1*b2)*g2*psi(idx,jdy)^2*conj(psi2(idx,jdy)))/(a2^2*b1^2)-(a1*sqrt(a2*b1)*b2*sqrt(a1*b2)*psi2(idx,jdy)^2*conj(psi2(idx,jdy)))/(a2^2*b1^2)+(2*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*psi(idx,jdy)*(ay(idx,jdy+1)-ay(idx,jdy))/dy)/(a2^2*b1*k1)-(1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*psi2(idx,jdy)*(ay(idx,jdy+1)-ay(idx,jdy))/dy)/(a2^2*b1*k1*m2)+(4*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*ay(idx,jdy)*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(a2^2*b1*k1)-(2*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*ay(idx,jdy)*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(a2^2*b1*k1*m2)-(2*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*(psi(idx,jdy+1)+psi(idx,jdy-1)-2*psi(idx,jdy))/dy^2)/(a2^2*b1*k1^2)+(a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*(psi2(idx,jdy+1)+psi2(idx,jdy-1)-2*psi2(idx,jdy))/dy^2)/(a2^2*b1*k1^2*m2)-(2*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*psi(idx,jdy)*(ax(idx+1,jdy)-ax(idx,jdy))/dx)/(a2^2*b1*k1)-(1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*psi2(idx,jdy)*(ax(idx+1,jdy)-ax(idx,jdy))/dx)/(a2^2*b1*k1*m2)-(4*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*ax(idx,jdy)*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(a2^2*b1*k1)-(2*1i*a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*ax(idx,jdy)*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(a2^2*b1*k1*m2)+(2*a1*sqrt(a2*b1)*sqrt(a1*b2)*g3*m1*(psi(idx+1,jdy)+psi(idx-1,jdy)-2*psi(idx,jdy))/dx^2)/(a2^2*b1*k1^2)+(a1*sqrt(a2*b1)*sqrt(a1*b2)*m1*(psi2(idx+1,jdy)+psi2(idx-1,jdy)-2*psi2(idx,jdy))/dx^2)/(a2^2*b1*k1^2*m2))*dt;
            psi2k(idx,jdy) = psi2(idx,jdy) + dpsi2;
            
            
            
            dpsi3 = (2*psi3(idx,jdy)-2*psi3(idx,jdy)^3-(a1*sqrt(a3/b3)*b3*l1*psi2(idx,jdy)*conj(psi(idx,jdy)))/(a3^2*b1)-(1i*a1*sqrt(a3/b3)*b3*l2*psi2(idx,jdy)*conj(psi(idx,jdy)))/(a3^2*b1)-(2*a1*l3*psi(idx,jdy)*psi3(idx,jdy)*conj(psi(idx,jdy)))/(a3*b1)-(a1*sqrt(a3/b3)*b3*l1*psi(idx,jdy)*conj(psi2(idx,jdy)))/(a3^2*b1)+(1i*a1*sqrt(a3/b3)*b3*l2*psi(idx,jdy)*conj(psi2(idx,jdy)))/(a3^2*b1)-(2*a1*l4*psi2(idx,jdy)*psi3(idx,jdy)*conj(psi2(idx,jdy)))/(a3*b1)+(2*a1*m1*(psi3(idx,jdy+1)+psi3(idx,jdy-1)-2*psi3(idx,jdy))/dy^2)/(a3*k1^2*m3)+(2*a1*m1*(psi3(idx+1,jdy)+psi3(idx-1,jdy)-2*psi3(idx,jdy))/dx^2)/(a3*k1^2*m3))*dt;
            psi3k(idx,jdy) = psi3(idx,jdy) + dpsi3;
                
                
           
            dax = (-((ax(idx,jdy)*psi(idx,jdy)*conj(psi(idx,jdy)))/ss)-(2*g3*m1*ax(idx,jdy)*psi2(idx,jdy)*conj(psi(idx,jdy)))/ss-(2*g3*m1*ax(idx,jdy)*psi(idx,jdy)*conj(psi2(idx,jdy)))/ss-(m1*ax(idx,jdy)*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(m2*ss)+(ax(idx,jdy+1)+ax(idx,jdy-1)-2*ax(idx,jdy))/dy^2/ss-(1i*conj(psi(idx,jdy))*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(2*k1*ss)-(1i*g3*m1*conj(psi2(idx,jdy))*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(k1*ss)-(1i*g3*m1*conj(psi(idx,jdy))*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(k1*ss)-(1i*m1*conj(psi2(idx,jdy))*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(2*k1*m2*ss)+(1i*psi(idx,jdy)*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(2*k1*ss)+(1i*g3*m1*psi2(idx,jdy)*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(k1*ss)+(1i*g3*m1*psi(idx,jdy)*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(k1*ss)+(1i*m1*psi2(idx,jdy)*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(2*k1*m2*ss)-(ay(idx+1,jdy+1)+ay(idx-1,jdy-1)-ay(idx-1,jdy+1)-ay(idx+1,jdy-1))/dx/dy/4/ss)*dt;
            axk(idx,jdy) = ax(idx,jdy) + dax;
            
            
            day = (-((ay(idx,jdy)*psi(idx,jdy)*conj(psi(idx,jdy)))/ss)+(2*g3*m1*ay(idx,jdy)*psi2(idx,jdy)*conj(psi(idx,jdy)))/ss+(2*g3*m1*ay(idx,jdy)*psi(idx,jdy)*conj(psi2(idx,jdy)))/ss-(m1*ay(idx,jdy)*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(m2*ss)-(1i*conj(psi(idx,jdy))*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(2*k1*ss)+(1i*g3*m1*conj(psi2(idx,jdy))*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(k1*ss)+(1i*g3*m1*conj(psi(idx,jdy))*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(k1*ss)-(1i*m1*conj(psi2(idx,jdy))*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(2*k1*m2*ss)+(1i*psi(idx,jdy)*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(2*k1*ss)-(1i*g3*m1*psi2(idx,jdy)*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(k1*ss)-(1i*g3*m1*psi(idx,jdy)*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(k1*ss)+(1i*m1*psi2(idx,jdy)*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(2*k1*m2*ss)-(ax(idx+1,jdy+1)+ax(idx-1,jdy-1)-ax(idx-1,jdy+1)-ax(idx+1,jdy-1))/dx/dy/4/ss+(ay(idx+1,jdy)+ay(idx-1,jdy)-2*ay(idx,jdy))/dx^2/ss)*dt;
            ayk(idx,jdy) = ay(idx,jdy) + day;
            
            err=err + abs(dpsi)^2;
            
            fre = fre -((a3^2*psi3(idx,jdy)^2)/b3)+(a3^2*psi3(idx,jdy)^4)/(2*b3)-(a1^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1+(a1^2*ax(idx,jdy)^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1+(a1^2*ay(idx,jdy)^2*psi(idx,jdy)*conj(psi(idx,jdy)))/b1+(2*a1^2*g3*m1*ax(idx,jdy)^2*psi2(idx,jdy)*conj(psi(idx,jdy)))/b1-(2*a1^2*g3*m1*ay(idx,jdy)^2*psi2(idx,jdy)*conj(psi(idx,jdy)))/b1+(a1*sqrt(a3/b3)*l1*psi2(idx,jdy)*psi3(idx,jdy)*conj(psi(idx,jdy)))/b1+(1i*a1*sqrt(a3/b3)*l2*psi2(idx,jdy)*psi3(idx,jdy)*conj(psi(idx,jdy)))/b1+(a1*a3*l3*psi(idx,jdy)*psi3(idx,jdy)^2*conj(psi(idx,jdy)))/(b1*b3)+(a1^2*psi(idx,jdy)^2*conj(psi(idx,jdy))^2)/(2*b1)+(a1^2*g2*psi2(idx,jdy)^2*conj(psi(idx,jdy))^2)/b1^2+(2*a1^2*g3*m1*ax(idx,jdy)^2*psi(idx,jdy)*conj(psi2(idx,jdy)))/b1-(2*a1^2*g3*m1*ay(idx,jdy)^2*psi(idx,jdy)*conj(psi2(idx,jdy)))/b1-(a1*a2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/b1+(a1^2*m1*ax(idx,jdy)^2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(b1*m2)+(a1^2*m1*ay(idx,jdy)^2*psi2(idx,jdy)*conj(psi2(idx,jdy)))/(b1*m2)+(a1*sqrt(a3/b3)*l1*psi(idx,jdy)*psi3(idx,jdy)*conj(psi2(idx,jdy)))/b1-(1i*a1*sqrt(a3/b3)*l2*psi(idx,jdy)*psi3(idx,jdy)*conj(psi2(idx,jdy)))/b1+(a1*a3*l4*psi2(idx,jdy)*psi3(idx,jdy)^2*conj(psi2(idx,jdy)))/(b1*b3)+(a1^2*g1*psi(idx,jdy)*psi2(idx,jdy)*conj(psi(idx,jdy))*conj(psi2(idx,jdy)))/b1^2+(a1^2*g2*psi(idx,jdy)^2*conj(psi2(idx,jdy))^2)/b1^2+(a1^2*b2*psi2(idx,jdy)^2*conj(psi2(idx,jdy))^2)/(2*b1^2)+(a1^2*(ax(idx,jdy+1)-ax(idx,jdy))/dy^2)/b1+(1i*a1^2*ay(idx,jdy)*conj(psi(idx,jdy))*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(b1*k1)-(2*1i*a1^2*g3*m1*ay(idx,jdy)*conj(psi2(idx,jdy))*(psi(idx,jdy+1)-psi(idx,jdy))/dy)/(b1*k1)-(2*1i*a1^2*g3*m1*ay(idx,jdy)*conj(psi(idx,jdy))*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(b1*k1)+(1i*a1^2*m1*ay(idx,jdy)*conj(psi2(idx,jdy))*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy)/(b1*k1*m2)+(a1*a3*m1*(psi3(idx,jdy+1)-psi3(idx,jdy))/dy^2)/(b3*k1^2*m3)-(1i*a1^2*ay(idx,jdy)*psi(idx,jdy)*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1)+(2*1i*a1^2*g3*m1*ay(idx,jdy)*psi2(idx,jdy)*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1)+(a1^2*(psi(idx,jdy+1)-psi(idx,jdy))/dy*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1^2)-(2*a1^2*g3*m1*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy*(conj(psi(idx,jdy+1))-conj(psi(idx,jdy)))/dy)/(b1*k1^2)+(2*1i*a1^2*g3*m1*ay(idx,jdy)*psi(idx,jdy)*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b1*k1)-(1i*a1^2*m1*ay(idx,jdy)*psi2(idx,jdy)*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b1*k1*m2)-(2*a1^2*g3*m1*(psi(idx,jdy+1)-psi(idx,jdy))/dy*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b1*k1^2)+(a1^2*m1*(psi2(idx,jdy+1)-psi2(idx,jdy))/dy*(conj(psi2(idx,jdy+1))-conj(psi2(idx,jdy)))/dy)/(b1*k1^2*m2)-(2*a1^2*(ax(idx,jdy+1)-ax(idx,jdy))/dy*(ay(idx+1,jdy)-ay(idx,jdy))/dx)/b1+(a1^2*(ay(idx+1,jdy)-ay(idx,jdy))/dx^2)/b1+(1i*a1^2*ax(idx,jdy)*conj(psi(idx,jdy))*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(b1*k1)+(2*1i*a1^2*g3*m1*ax(idx,jdy)*conj(psi2(idx,jdy))*(psi(idx+1,jdy)-psi(idx,jdy))/dx)/(b1*k1)+(2*1i*a1^2*g3*m1*ax(idx,jdy)*conj(psi(idx,jdy))*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(b1*k1)+(1i*a1^2*m1*ax(idx,jdy)*conj(psi2(idx,jdy))*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx)/(b1*k1*m2)+(a1*a3*m1*(psi3(idx+1,jdy)-psi3(idx,jdy))/dx^2)/(b3*k1^2*m3)-(1i*a1^2*ax(idx,jdy)*psi(idx,jdy)*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1)-(2*1i*a1^2*g3*m1*ax(idx,jdy)*psi2(idx,jdy)*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1)+(a1^2*(psi(idx+1,jdy)-psi(idx,jdy))/dx*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1^2)+(2*a1^2*g3*m1*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx*(conj(psi(idx+1,jdy))-conj(psi(idx,jdy)))/dx)/(b1*k1^2)-(2*1i*a1^2*g3*m1*ax(idx,jdy)*psi(idx,jdy)*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b1*k1)-(1i*a1^2*m1*ax(idx,jdy)*psi2(idx,jdy)*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b1*k1*m2)+(2*a1^2*g3*m1*(psi(idx+1,jdy)-psi(idx,jdy))/dx*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b1*k1^2)+(a1^2*m1*(psi2(idx+1,jdy)-psi2(idx,jdy))/dx*(conj(psi2(idx+1,jdy))-conj(psi2(idx,jdy)))/dx)/(b1*k1^2*m2);
            
        end
    end
    
	EEE2=[EEE2,fre/Nx/Ny];
    
   


    
    
    

    err=err/Nx/Ny;
    dt=0.005*(1-tanh((err - 1e-7)/0.01));

    if err<errk
        errk=err;
    elseif count>1000
        fre=0;
        break;
    end
    
    errk=err;
    
    ax = axk;
    ay = ayk;
    psi = psik;
    psi2 = psi2k;
    psi3 = psi3k;
    count=count+1;
    if mod(count,ddt)==0
       err
        m=matfile(sprintf('MYGL%d%d.mat', lda,count),'writable',true);
        m.x=psi;
        m.y=psi2;
        m.z=psi3;
        m.l=ax;
        m.n=ay;
    end

end

for jdy = 2:Ny+1
        psi(Nx+2,jdy) = psi(2,jdy)*exp(+1i*(jdy-1)*dy*phi/2/Ly);
        psi(1,jdy) = psi(Nx+1,jdy)*exp(-1i*(jdy-1)*dy*phi/2/Ly);
        psi2(Nx+2,jdy) = psi2(2,jdy)*exp(+1i*(jdy-1)*dy*phi/2/Ly);
        psi2(1,jdy) = psi2(Nx+1,jdy)*exp(-1i*(jdy-1)*dy*phi/2/Ly);
        psi3(Nx+2,jdy) = psi3(2,jdy);
        psi3(1,jdy) = psi3(Nx+1,jdy);
        
        ax(Nx+2,jdy) = ax(2,jdy);
        ax(1,jdy) = ax(Nx+1,jdy);
        ay(Nx+2,jdy) = ay(2,jdy) + phi/2/Ly/k1;
        ay(1,jdy) = ay(Nx+1,jdy) - phi/2/Ly/k1;
end

    for idx = 2:Nx+1
        psi(idx,Ny+2) = psi(idx,2)*exp(-1i*(idx-1)*dx*phi/2/Lx);
        psi(idx,1) = psi(idx,Ny+1)*exp(+1i*(idx-1)*dx*phi/2/Lx);
        psi2(idx,Ny+2) = psi2(idx,2)*exp(-1i*(idx-1)*dx*phi/2/Lx);
        psi2(idx,1) = psi2(idx,Ny+1)*exp(+1i*(idx-1)*dx*phi/2/Lx);
        psi3(idx,Ny+2) = psi3(idx,2);
        psi3(idx,1) = psi3(idx,Ny+1);
        
        ax(idx,Ny+2) = ax(idx,2) - phi/2/Lx/k1;
        ax(idx,1) = ax(idx,Ny+1) + phi/2/Lx/k1;
        ay(idx,Ny+2) = ay(idx,2);
        ay(idx,1) = ay(idx,Ny+1);
    end
    
    ax(Nx+2,Ny+2) = ax(2,2)-phi/2/Lx/k1;
    ax(1,1) = ax(Nx+1,Ny+1)+phi/2/Lx/k1;
    ay(Nx+2,Ny+2) = ay(2,2)+phi/2/Ly/k1;
    ay(1,1) = ay(Nx+1,Ny+1)-phi/2/Ly/k1;


EEE3=[EEE3,fre/Nx/Ny];

        m=matfile(sprintf('MYGLapppttb%d.mat', lda),'writable',true);
        m.x=psi;
        m.y=psi2;
        m.z=psi3;
        m.l=ax;
        m.n=ay;
        
        me=matfile(sprintf('MYGLeeeb2.mat'),'writable',true);
        me.u=EEE2;
        me.v=EEE3;
    







end

%%% Plot this
% Discretize the X and Y space
% x_vec = linspace(0,Lx,Nx);
% y_vec = linspace(0,Ly,Ny);

fix(clock)
xxx=1:length(EEE3);
plot(xxx,real(EEE3))

% psip = zeros(Nx,Ny);
% for idx = 1:Nx
%     for jdy = 1:Ny
%         psip(idx,jdy)=pp(idx+1,jdy+1,floor(Nt/ddt))*conj(pp(idx+1,jdy+1,floor(Nt/ddt)));
%     end
% end 
% 
% 
% 
% [yy,xx] = meshgrid(y_vec,x_vec);
% contourf(yy,xx,psip)
% set(gca,'DataAspectRatio',[1 1 1])
% xlabel('X/m')
% ylabel('Y/m')
% zlabel('Psi')




