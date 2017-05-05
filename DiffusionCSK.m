% 2D Transient Heat Equation solver via finite-difference scheme

clear;
clc;
Lx=2*pi;
Ly=2*pi;
nx=41;
ny=41;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
dt=.01; 
t_f=2;
T0=1;
T1= 0;
% Initial condition
r1=0.5*dt/(dx^2);
r2=0.5*dt/(dy^2);

if (r1>=0.5)
    error('Unstable Solution!Please change dx dy and dt.');
end
if (r2>=0.5)
    error('Unstable Solution!Please change dx dy and dt.');
end


for i=1:nx
    for j=1:ny
        T(i,j,1)=T0;
    end
end

% Boundary conditions
ax=0;ay=0;bx=2*pi;by=2*pi;
for i=1:nx
    T(i,1,1)=g(2*pi*(i-1)/nx);
    T(i,ny,1)=f(2*pi*(i-1)/nx);
end
alpha=0.5;
for j=1:ny
    T(1,j,1)=g(bx)-(2*pi*(j-1)/ny-ay)*(f(bx)-g(bx))/(by-ay);
    T(nx,j,1)=0;
end

for t=1:((t_f/dt))
    
    for i=1:nx
        T(i,1,t)=0;
        T(i,ny,t)=0;
    end

    for j=1:ny
        T(1,j,t)=0;
        T(nx,j,t)=0;
    end
    
    for it=1:100
        for i=2:(nx-1)
            for j=2:(ny-1)
                T(i,j,t+1)=alpha*dt*((T(i+1,j,t)-2*T(i,j,t)+T(i-1,j,t))/(dx^2)+(T(i,j+1,t)-2*T(i,j,t)+T(i,j-1,t))/(dy^2))+T(i,j,t);
            end
        end
    end
    x=linspace(0,Lx,nx);
    y=linspace(0,Ly,ny);
    surf(y,x,T(:,:,t+1))
    title(['Heat Transfer at t(seconds) = ',num2str((t*dt))])
    axis([0 Ly 0 Lx 0 1]);
    eval(['print -djpeg heat2d_' num2str(t) '.jpeg']);
end

