clear, close;

%% defining the constants
lambda    = 1;
epsilonr1  = 4;
epsilonr2  = 4;
k0        = 2*pi/lambda;
k1        = k0*sqrt(epsilonr1);
k2        = k0*sqrt(epsilonr2);
phi0      = 0; %incident angle
theta0    = pi/2;
E0        = 1; %incident field amplitude
gamma     = 0.577; %for hankel fun approximation
n         = 2; %for gaussquad

%% building the scatterers
%here we will make two cylindrical discs in the domain
N          = 50; %number of segments on the surface(boundary)
radius     = [0.3*lambda, 0.4*lambda];
circum     = 2*pi.*radius;
l          = circum/N; %length of each segment
step       = 2*pi/N; %angle measure of each segment
nodes      = 0:step:2*pi-step; %nodes of each segment
%testing points which are choosen as the center of each segment
test_pts   = step/2 : step : 2*pi-step/2; 

%scatterer 1
scatter(-lambda  + radius(1)*cos(nodes),radius(1)*sin(nodes),'k','filled');
hold on; grid on; axis('equal');
%scatterer 2
scatter(lambda + radius(2)*cos(nodes),radius(2)*sin(nodes),'k','filled');
% pause(2)
% scatter(-lambda + radius(1)*cos(test_pts),radius(1)*sin(test_pts),'r','filled');
% scatter(lambda + radius(2)*cos(test_pts),radius(2)*sin(test_pts),'r','filled');
% grid on; axis('equal');

%% formulating the problem
%defining the known side of Ax = c i.e. c
%for phi_i(or Ei) vector

%scatterer 1 coordinates
X1 = -lambda + radius(1)*cos(test_pts);
Y1 = radius(1)*sin(test_pts);
%scatterer 2 coordinates
X2 = lambda + radius(2)*cos(test_pts);
Y2 = radius(2)*sin(test_pts);

%angle for the incident plane wave
%for scatterer 1
alpha1 = (X1 * sin(theta0) * cos(phi0)) + (Y1 * sin(theta0) * sin(phi0));
%for scatterer 2
alpha2 = (X2 * sin(theta0) * cos(phi0)) + (Y2 * sin(theta0) * sin(phi0));

%incident electric field on surface of scatterer 1
Ei1 = E0 * exp(-1i * k1 * alpha1);     
%incident electric field on surface of scatterer 2
Ei2 = E0 * exp(-1i * k2 * alpha2);     


%% plotting the result

%% defining functions
%% function definitions
%%green function
function g = green(rpnl,rpnu,n,rnl,rnu,m,k)
    rps = ((1-n).*rpnl) + (n.*rpnu);
    rt  = ((1-m).*rnl)  + (m.*rnu);
    rho = norm(rt-rps);
    g   = (-1j/4) * besselh(0,2,k.*rho);
end

%%grad_green
function gg = gradg(rpnl,rpnu,n,rnl,rnu,m,nhat,k)
    rps  = ((1-n).*rpnl) + (n.*rpnu);
    rt   = ((1-m).*rnl)  + (m.*rnu);
    rho  = norm(rt-rps);
    rhat = (rt-rps)./rho;
    gg   = (1j*k/4) * besselh(1,2,k.*rho) .* dot(rhat,nhat);
end

%%gauss quadrature
function y = gaussquad(f,n)
%     wt = zeros(n,1);
%     xi = zeros(n,1);
    if n == 1
        wt(1)  = 1; 
        xi(1) = 1;
    elseif n == 2
        xi(1) = (1 - 1/sqrt(3))/2;
        xi(2) = (1 + 1/sqrt(3))/2;
        wt(1) = 0.5;
        wt(2) = 0.5;
    end
    y = 0;
    for i = 1:n
        y = y + wt(i) * f(xi(i));
    end
end