% Solve 1D Euler equation
% n - number of spatial grids
% cfl - CFL number
% t_F - final time
% rhol - initial density on the left
% rhor - initial density on the right
% rhoul - initial momentum on the left
% rhour - initial momentum on the right
% rhoel - initial total energy on the left
% rhoer - initial total energy on the right
% called by MATLAB command line: kfvs_Euler
% Author: Haohan Li, hlibb@connect.ust.hk
function kfvs_Euler

%% initial condition
n = 500; cfl = 0.4; t_F = 0.2; rhol = 1.0; rhor = 0.1;
rhoul = 0.0; rhour = 0.0; rhoel = 1.0; rhoer = 0.25;
xl = 0.0; xr = 1.0; m = 1; 
gammar = 1.4;
dx = (xr-xl)/n;
x = zeros(n);
w0 = zeros(n+2*m,3);
w1 = zeros(n+2*m,3);
rho = zeros(n);
u = zeros(n);
p = zeros(n);
c = zeros(n);
[x_exact,rho_exact,u_exact,p_exact,e_exact]=textread('./Euler1D/e1rpex.out','%f%f%f%f%f','headerlines',1);

for i = 1:n
    x(i) = xl+0.5*dx+(i-1)*dx;
end
 
for i = 1:n
    if x(i) <= 0.5*(xr+xl)
        w0(i+m,1) = rhol;
        w0(i+m,2) = rhoul;
        w0(i+m,3) = rhoel;
    else
        w0(i+m,1) = rhor;
        w0(i+m,2) = rhour;
        w0(i+m,3) = rhoer;
    end
end

%% boundary condition
for i = 1:m
    w0(i,:) = w0(1+m,:);
    w0(n+m+i,:) = w0(n+m,:);
end

%% calculate density, velocity and pressure
for i = 1:n
	rho(i) = w0(m+i,1);
	u(i) = w0(m+i,2)/w0(m+i,1);	
	p(i) = (gammar-1.)*(w0(m+i,3)-0.5*rho(i)*u(i)*u(i));
end

%% plot initial
%% plot initial
figure
hold on
plot(x,rho,'b--','LineWidth',1); 
plot(x_exact,rho_exact,'b-','LineWidth',2);  
xlabel('x'); ylabel('density');
axis([0,1,0.0,1.1]);
hold off

figure
hold on
plot(x,u,'k--','LineWidth',1); 
plot(x_exact,u_exact,'k-','LineWidth',2);  
xlabel('x'); ylabel('velocity');
axis([0,1,-0.05,0.55]);
hold off

figure
hold on
plot(x,p,'r--','LineWidth',1); 
plot(x_exact,p_exact,'r-','LineWidth',2); 
xlabel('x'); ylabel('pressure');
axis([0,1,0.0,0.45]);
hold off
pause;

%% start iteration
cnt = 0;
t = 0.0;
while cnt <= 1000000 % maximum steps
    cnt = cnt+1;
    sigmamax = 0.0;
    for i = 1:n
        rho0 = w0(m+i,1);
        u0 = w0(m+i,2)/w0(m+i,1);	
        p0 = (gammar-1.)*(w0(m+i,3)-0.5*rho(i)*u(i)*u(i));
        c = sqrt(gammar*p0/rho0); % sound speed
        sigma = max(abs(u0+c),abs(u0-c));
	    sigmamax = max(sigma,sigmamax);	
    end
    dt = cfl*dx/sigmamax; % CFL condition
    if t < t_F && t+dt > t_F
        dt = t_F-t;
    end
    t = t+dt % time stepping
    if t >= t_F
        for i = 1:n
            rho(i) = w0(m+i,1);
	        u(i) = w0(m+i,2)/w0(m+i,1);	
	        p(i) = (gammar-1.)*(w0(m+i,3)-0.5*rho(i)*u(i)*u(i));
        end
        % plot solution
        %% plot initial
        figure
        hold on
        plot(x,rho,'b--','LineWidth',1); 
        plot(x_exact,rho_exact,'b-','LineWidth',2);  
        xlabel('x'); ylabel('density');
        axis([0,1,0.0,1.1]);
        hold off

        figure
        hold on
        plot(x,u,'k--','LineWidth',1); 
        plot(x_exact,u_exact,'k-','LineWidth',2);  
        xlabel('x'); ylabel('velocity');
        axis([0,1,-0.05,0.55]);
        hold off

        figure
        hold on
        plot(x,p,'r--','LineWidth',1); 
        plot(x_exact,p_exact,'r-','LineWidth',2); 
        xlabel('x'); ylabel('pressure');
        axis([0,1,0.0,0.45]);
        hold off
        break;
    end
%% update W
    for i = 1:n
        w1(m+i,:) = w0(m+i,:)+dt/dx...
      *(k_f_Euler(w0(m+i-1,:),w0(m+i,:))-k_f_Euler(w0(m+i,:),w0(m+i+1,:)));
    end
    % boundary condition
    for i = 1:m
        w1(i,:) = w1(1+m,:);
        w1(n+m+i,:) = w1(n+m,:);
    end
    % update
    w0(:,:,:) = w1(:,:,:);
    % every 100 steps draw results out
    if mod(cnt,300) == 0
        for i = 1:n
            rho(i) = w0(m+i,1);
	        u(i) = w0(m+i,2)/w0(m+i,1);	
	        p(i) = (gammar-1.)*(w0(m+i,3)-0.5*rho(i)*u(i)*u(i));
        end
        % plot solution
        %% plot initial
        figure
        hold on
        plot(x,rho,'b--','LineWidth',1); 
        plot(x_exact,rho_exact,'b-','LineWidth',2);  
        xlabel('x'); ylabel('density');
        axis([0,1,0.0,1.1]);
        hold off

        figure
        hold on
        plot(x,u,'k--','LineWidth',1); 
        plot(x_exact,u_exact,'k-','LineWidth',2);  
        xlabel('x'); ylabel('velocity');
        axis([0,1,-0.05,0.55]);
        hold off

        figure
        hold on
        plot(x,p,'r--','LineWidth',1); 
        plot(x_exact,p_exact,'r-','LineWidth',2); 
        xlabel('x'); ylabel('pressure');
        axis([0,1,0.0,0.45]);
        hold off
        pause;
    end
end


    
