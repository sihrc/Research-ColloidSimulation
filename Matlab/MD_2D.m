%% Molecular Dynamics Simulation
% <MD_2D.m> Thibault Bertrand 2/23/2013

% Based on <Vibr2DCell.m> by Thibault Bertrand
% MD operated by Velocity-Verlet
% Forces: 
% ======
%       -interparticle interaction, linear elastic repulsion
%       -wall interaction, linear elastic repulsion
%       -collisional damping
%       -brownian thermostat
%       -Fluid damping
%       
% Revision history:
% 11/2/2012: added wall interaction
% 11/2/2012: added polydispersity via polydispersity parameter pp
% 11/2/2012: added Kw, spring constant for interaction with the walls
% 2/23/2013: added the fluid damping and the brownian thermostat through RNG 


%% Experimental Parmeters
tic
Ns=200;  % number of small particles
Nl=0;  % number of large particles (if any bidispersity)
D=2;    % small diameter of particles
G=1.4;  % ratio of large to small diameter
K=100;  % spring constant for harmonic force law 
Kw=500; % spring constant for interaction with the walls 
M=3;    % mass of particles
gg=1;   % collisional damping
pp=0.05; % polydispersity parameter

B=0.1;  % air drag damping
Lx=50*D; % width of box
Ly=20*D; % height of box

%% Physical Parameters
g=0.1;

%% calculted parameters
N=Ns+Nl;
Gn=[ones(1,Ns) G*ones(1,Nl)].*(1+pp*randn(1,Ns+Nl)); % w polydispersity, polydispersity paramter pp
%Gn=[ones(1,Ns) G*ones(1,Nl)]; % w/o polydispersity
Dn=D*Gn;
xlp=-Lx/2; % used in no PBC in x
xrp=Lx/2; % used in no PBC in x
ybp=0;
ytp=Ly;

%% Display Parameters
Nplotskip=20;  % number of timesteps to skip before plotting

%% Simulation Parmeters
dt=1e-2; % time step of the simulation
Ne=30000; % Number of steps (equilibration)
Nt=10000; % Number of steps (recorded simulation)

%% Random Generator Number Parameters
RNG_on=1;
T=1e-1; % temperature (controls the RMSD of the particles through the RNG)
eta=1; % viscosity of the fluid
Diff=T/(6*pi*eta*D/2);
sigma=sqrt(2*Diff*dt); % standard deviation of the gaussian distribution for RNG
NRNGskip=100;
rng('shuffle'); % Seeds the RNG with the clock time as a seed so you get a different sequence each time

%% Saved Parameters
Ek=zeros(1,Nt+Ne);

%% Initial Conditions
[x y]=ndgrid(-Lx/2+G*D/2:G*D:Lx/2-G*D/2,G*D/2:G*D:Ly-G*D/2);  
[junk ii]=sort(rand(1,numel(x)));

x=x(ii(1:N));
y=y(ii(1:N));
vx=randn(1,N)/1000;
vy=randn(1,N)/1000;
ax_old=0*x;
ay_old=0*y;

%% Setup Plotting
% clf;
% hb=plot([-Lx/2,Lx/2],[0,0]);
% hold on;
% ht=plot([-Lx/2,Lx/2],[Ly,Ly]);
% h=zeros(1,N);
% for np=1:N
%   h(np)=rectangle('Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
% end
% axis('equal');
% axis([-Lx/2 Lx/2 0 Ly]);
% hold off;
%pause;


%% Main Loop

for nt=1:Nt+Ne
  
  % Generation of random numbers
  if(RNG_on == 1 && rem(nt-1,NRNGskip)==0)
     xrand=sigma*randn(N,NRNGskip); % Random numbers on a gaussian distribution 
     yrand=sigma*randn(N,NRNGskip); % of mean 0 and standard deviation 
  end
  
  x=x+xrand(:,rem(nt-1,NRNGskip)+1)';
  y=y+yrand(:,rem(nt-1,NRNGskip)+1)';
  
%   % plot particles
%   if(rem(nt-1,Nplotskip)==0)
%     for np=1:N
%       set(h(np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
%     end
%     if(nt<Ne)
%         title(['N=' num2str(N) '  T=' num2str(T) '-- Equilibration']);
%     else
%         title(['N=' num2str(N) '  T=' num2str(T) '-- Steady State']);    
%     end
%     drawnow;
%   end
%   
  % First step in Velocity Verlet integration
  x=x+vx*dt+ax_old.*dt.^2/2;  
  y=y+vy*dt+ay_old.*dt.^2/2;
  
  Ek(nt)=sum(M.*(vx.^2+vy.^2))/2;
  
  %x=x-Lx*round(x/Lx);  % Periodic x -- comment out if not PBC in x
  
  % Interaction detector and Force Law
  Fx=zeros(1,N);
  Fy=zeros(1,N);
  
  for nn=1:N
    for mm=nn+1:N
      dy=y(mm)-y(nn);
      Dnm=(Dn(mm)+Dn(nn))/2;
      if(dy<Dnm)
        dx=x(mm)-x(nn);
        %dx=dx-round(dx/Lx)*Lx;  % Periodic x -- comment out if not PBC in x
        dnm=dx.^2+dy.^2;
        if(dnm<Dnm^2)
          dnm=sqrt(dnm);
          F=-K*(Dnm/dnm-1);
          dvx=vx(mm)-vx(nn);
          dvy=vy(mm)-vy(nn);
          Fx(nn)=Fx(nn)+F.*dx+gg*dvx;  % particle-particle Force Law
          Fx(mm)=Fx(mm)-F.*dx-gg*dvx;
          Fy(nn)=Fy(nn)+F.*dy+gg*dvy;  % particle-particle Force Law
          Fy(mm)=Fy(mm)-F.*dy-gg*dvy;
        end
      end
    end
  end
  
  ii=(y<ybp+Dn/2);
  Fy(ii)=Fy(ii)-Kw*(y(ii)-Dn(ii)/2-ybp);  % Bottom wall

  ii=(y>ytp-Dn/2);
  Fy(ii)=Fy(ii)-Kw*(y(ii)-(ytp-Dn(ii)/2));  % top wall
  
  ii=(x<xlp+Dn/2);
  Fx(ii)=Fx(ii)-Kw*(x(ii)-Dn(ii)/2-xlp);  % Left wall --  uncommented if not PBC in x
  
  ii=(x>xrp-Dn/2);
  Fx(ii)=Fx(ii)-Kw*(x(ii)+Dn(ii)/2-xrp);  % Right wall --  uncommented if not PBC in x

 
  ax=(Fx-B*vx)./M;
  ay=(Fy-B*vy)./M-g;
  
  % Second step in Verlet integration
  vx=vx+(ax_old+ax).*dt/2;  
  vy=vy+(ay_old+ay).*dt/2;

  ax_old=ax;
  ay_old=ay;
  
end
toc