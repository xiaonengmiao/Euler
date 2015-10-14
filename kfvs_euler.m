% Solve 1D Euler equation
% nx - number of spatial grids
% cfl - CFL number
% t_F - final time
% rhol - initial density on the left
% rhor - initial density on the right
% rhoul - initial momentum on the left
% rhour - initial momentum on the right
% rhoEl - initial total energy on the left
% rhoEr - initial total energy on the right
% rk -  1, Kinetic flux  2, Roe flux 
% fre - output frequency
% called by MATLAB command line:    kfvs_euler(500, 0.4, 0.2, 1., 0.1, 0., 0., 1., 0.25, 1, 100)
% Author: Luo Li, lluoac@ust.hk
function kfvs_euler(nx,cfl,t_F,rhol,rhor,rhoul,rhour,rhoEl,rhoEr,rk,fre)
gammer=1.4;
itmax=100000; % maximum of time steps
m=1; % number of ghost points in x
xmin=0.;
xmax=1.;
dx=(xmax-xmin)/nx; % space step length
x=zeros(nx);
u0=zeros(nx+2*m,3);  
u1=zeros(nx+2*m,3);
rho_=zeros(nx);  
u_=zeros(nx); 
p_=zeros(nx); 

[x_exact,rho_exact,u_exact,p_exact,e_exact]=textread('./Euler1D/e1rpex.out','%f%f%f%f%f','headerlines',1);

for i=1:nx
    x(i)=xmin+0.5*dx+(i-1)*dx;
end

% initialize
    for i=1:nx
        if x(i)<=0.5*(xmax-xmin)
           u0(m+i,1)=rhol;
           u0(m+i,2)=rhoul;
           u0(m+i,3)=rhoEl;		   
        else
           u0(m+i,1)=rhor;
           u0(m+i,2)=rhour;
           u0(m+i,3)=rhoEr;	
        end
    end

% apply boundary condition
	for i=1:m
		% left outflow
		u0(i,:)=u0(m+1,:);
		% right outflow
		u0(m+nx+i,:)=u0(m+nx,:);
	end


% plot
for i=1:nx
	rho_(i)=u0(m+i,1);
	u_(i)=u0(m+i,2)/u0(m+i,1);	
	p_(i)=(gammer-1.)*(u0(m+i,3)-0.5*u0(m+i,2)*u0(m+i,2)/u0(m+i,1));
end
plot(x,rho_,'b--','LineWidth',1); hold on;
plot(x_exact,rho_exact,'b-','LineWidth',2); hold on; 
xlabel('x'); ylabel('density');
axis([0,1,-0.05,1.2]);
figure;
plot(x,u_,'k--','LineWidth',1); hold on;
plot(x_exact,u_exact,'k-','LineWidth',2); hold on; 
xlabel('x'); ylabel('velocity');
axis([0,1,-0.05,0.55]);
figure;
plot(x,p_,'r--','LineWidth',1); hold on;
plot(x_exact,p_exact,'r-','LineWidth',2); hold on;
xlabel('x'); ylabel('pressure');
axis([0,1,-0.05,0.45]);
pause;

% start time stepping
t=0.0;
for it=1:itmax
	rhomaxx=0.;
		for i=1:nx
			   % primitive variables
			   rho=u0(m+i,1);
			   rhou=u0(m+i,2);
			   rhoE=u0(m+i,3);
			   u=rhou/rho;
			   p=(gammer-1.)*(rhoE-0.5*rhou*rhou/rho);
			   % sound speed
			   c = sqrt(gammer*p/rho);
			   % evaluate the maximum of velocity
			   rhox = max(abs(u+c),abs(u-c));
			   rhomaxx = max(rhomaxx,rhox);		   
		end

    % compute time step length 
    dt=cfl*dx/rhomaxx;
	
    % final time step length
    if t<t_F && t+dt>t_F
        dt=t_F-t;
    end
    t=t+dt
	
    % if reach time limit, stop time stepping
    if t+dt>t_F || it==itmax
		% plot
		for i=1:nx
			rho_(i)=u0(m+i,1);
			u_(i)=u0(m+i,2)/u0(m+i,1);	
			p_(i)=(gammer-1.)*(u0(m+i,3)-0.5*u0(m+i,2)*u0(m+i,2)/u0(m+i,1));
		end
		plot(x,rho_,'b--','LineWidth',1); hold on;
		plot(x_exact,rho_exact,'b-','LineWidth',2); hold on; 
		xlabel('x'); ylabel('density');
		axis([0,1,-0.05,1.2]);
		figure;
		plot(x,u_,'k--','LineWidth',1); hold on;
		plot(x_exact,u_exact,'k-','LineWidth',2); hold on;
		xlabel('x'); ylabel('velocity');
        axis([0,1,-0.05,0.55]);
		figure;
		plot(x,p_,'r--','LineWidth',1); hold on;
		plot(x_exact,p_exact,'r-','LineWidth',2); hold on;
		xlabel('x'); ylabel('pressure');
		axis([0,1,-0.05,0.45]);
        break;
    end
    
    lambdax=dt/dx;
	
    % update solution
	if rk==1
		for i=1:nx
			u1(m+i,:)=u0(m+i,:)+...
				lambdax*(k_f_Euler(u0(m+i-1,:),u0(m+i,:))-k_f_Euler(u0(m+i,:),u0(m+i+1,:)));
		end
	else
		for i=1:nx
			u1(m+i,:)=u0(m+i,:)+...
				lambdax*(r_f_euler(u0(m+i-1,:),u0(m+i,:))-r_f_euler(u0(m+i,:),u0(m+i+1,:)));
		end	
	end
	
	% apply boundary condition
	for i=1:m
		% left outflow
		u1(i,:)=u1(m+1,:);
		% right outflow
		u1(m+nx+i,:)=u1(m+nx,:);
	end
	
    % swap
    u0(:,:,:)=u1(:,:,:);
	
    if mod(it,fre)==0 % draw frequency
		for i=1:nx
			rho_(i)=u0(m+i,1);
			u_(i)=u0(m+i,2)/u0(m+i,1);	
			p_(i)=(gammer-1.)*(u0(m+i,3)-0.5*u0(m+i,2)*u0(m+i,2)/u0(m+i,1));
		end
		plot(x,rho_,'b--','LineWidth',1); hold on;
		plot(x_exact,rho_exact,'b-','LineWidth',2); hold on; 
		xlabel('x'); ylabel('density');
		axis([0,1,-0.05,1.2]);
		figure;
		plot(x,u_,'k--','LineWidth',1); hold on;
		plot(x_exact,u_exact,'k-','LineWidth',2); hold on;
		xlabel('x'); ylabel('velocity');
        axis([0,1,-0.05,0.55]);
		figure;
		plot(x,p_,'r--','LineWidth',1); hold on;
		plot(x_exact,p_exact,'r-','LineWidth',2); hold on;
		xlabel('x'); ylabel('pressure');
		axis([0,1,-0.05,0.45]);
		pause;
    end
end
