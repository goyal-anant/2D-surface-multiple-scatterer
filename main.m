clear, close;
tic
%% defining the constants
lambda    = 1;
epsilonr1 = 4;
epsilonr2 = 4;
k0        = 2*pi/lambda;
k1        = k0*sqrt(epsilonr1);
k2        = k0*sqrt(epsilonr2);
phi0      = 0; %incident angle
theta0    = pi/2;
E0        = 1; %incident field amplitude
gamma     = 0.577; %for hankel fun approximation
n         = 4; %for gaussquad

%% building the scatterers
%here we will make two cylindrical discs in the domain
N          = 50; %number of segments on the surface(boundary)
radius     = [0.3*lambda, 0.3*lambda];
circum     = 2*pi.*radius;
l          = circum/N; %length of each segment
step       = 2*pi/N; %angle measure of each segment
nodes      = 0:step:2*pi-step; %nodes of each segment
%testing points which are choosen as the center of each segment
test_pts   = step/2 : step : 2*pi-step/2; 

%scatterer 1
% scatter(-lambda  + radius(1)*cos(nodes),radius(1)*sin(nodes),'k','filled');
% hold on; grid on; axis('equal');
%scatterer 2
% scatter(lambda + radius(2)*cos(nodes),radius(2)*sin(nodes),'k','filled');
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

%creating the A matrix
%A = [Q T U V; W X O O; O O Y Z]
Q = zeros(2*N,N); T = Q; U = Q; V = Q;
W = zeros(N,N);   X = W; Y = W; Z = W;

%integrable singularity, i.e., int(phi)
gdiag = @(k) -1j/4*(l(1) - k^2*l(1)^3/48 - 1j*(2*l(1)/pi*(log(l(1)*k/(4*exp(1)))+gamma)));

%loop for r prime
for i = 1:2*N
    if i <= N
        %first object
        rpnl = [-lambda + radius(1)*cos(test_pts(i) - (step/2)), radius(1)*sin(test_pts(i) - (step/2))];
        rpnu = [-lambda + radius(1)*cos(test_pts(i) + (step/2)), radius(1)*sin(test_pts(i) + (step/2))];
    else
        %second object
        rpnl = [lambda + radius(2)*cos(test_pts(i-N) - (step/2)), radius(2)*sin(test_pts(i-N) - (step/2))];
        rpnu = [lambda + radius(2)*cos(test_pts(i-N) + (step/2)), radius(2)*sin(test_pts(i-N) + (step/2))];
    end
    %loop for r
    for j = 1:2*N
        if j <= N
            %first object
            rnl = [-lambda + radius(1)*cos(test_pts(j) - (step/2)), radius(1)*sin(test_pts(j) - (step/2))];
            rnu = [-lambda + radius(1)*cos(test_pts(j) + (step/2)), radius(1)*sin(test_pts(j) + (step/2))];
            nhat= [-lambda + cos(test_pts(j)), sin(test_pts(j))];
        else
            %second object
            rnl = [lambda + radius(2)*cos(test_pts(j-N) - (step/2)), radius(2)*sin(test_pts(j-N) - (step/2))];
            rnu = [lambda + radius(2)*cos(test_pts(j-N) + (step/2)), radius(2)*sin(test_pts(j-N) + (step/2))];
            nhat= [lambda + cos(test_pts(j-N)), sin(test_pts(j-N))];
        end
        
        %green and grad green for region 0
        gk0  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k0);
        ggk0 = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k0);
        %green and grad green for region 1
        gk1  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k1);
        ggk1 = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k1);
        %green and grad green for region 2
        gk2  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k2);
        ggk2 = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k2);

        %r on object 1 and r prime on object 1
        if i<=N && j<=N
            if i==j
                Q(i,j) = gdiag(k0);
                T(i,j) = 1/2;
            else
                Q(i,j) =  l(1) * gaussquad(gk0,n);
                T(i,j) = -l(1) * gaussquad(ggk0,n);
            end
        %r on object 2 and r prime on object 1
        elseif i<=N && j>N
            U(i,j-N) = -l(2) * gaussquad(gk0,n);
            V(i,j-N) = -l(2) * gaussquad(ggk0,n);

            Y(i,j-N) =  l(2) * gaussquad(gk2,n);
            Z(i,j-N) = -l(2) * gaussquad(gk2,n);
        %r on object 1 and r prime on object 2
        elseif i>N && j<=N
            Q(i,j)   =  l(1) * gaussquad(gk0,n);
            T(i,j)   = -l(1) * gaussquad(ggk0,n);

            W(i-N,j) =  l(1) * gaussquad(gk1,n);
            X(i-N,j) = -l(1) * gaussquad(gk1,n);
        %r on object 2 and r prime on object 2
        elseif i>N && j>N
            if i==j
                %check if there would be -ve sign here or not
                U(i,j-N) = -gdiag(k0);
                V(i,j-N) = -1/2;
            else
                U(i,j-N) = -l(1) * gaussquad(gk0,n);
                V(i,j-N) = -l(1) * gaussquad(ggk0,n);
            end
        end
    end
end
%A = [Q T U V; W X O O; O O Y Z]
A = [Q, T, U, V; W, X, zeros(N,N), zeros(N,N); zeros(N,N), zeros(N,N), Y, Z];

%creating y in Ax = y
%y(1:N)=incident field on S1
%y(N+1:2N)=incident field on S2
y = [Ei1, Ei2, zeros(1,2*N)];


%solution vector
x = A\y';
% imagesc(abs(A)); colorbar;

%finding the total field now using huygen principle
oradius   = 200 * lambda;
oN        = 100;
ostep     = 2*pi/oN;
onodes    = 0:ostep:2*pi-ostep;
otest_pts = ostep/2 : ostep : 2*pi-ostep/2;
farfield  = zeros(1,oN); %scattered far field

for i = 1:oN
    rpnl = [oradius*cos(otest_pts(i) - (ostep/2)) oradius*sin(otest_pts(i) - (ostep/2))];
    rpnu = [oradius*cos(otest_pts(i) + (ostep/2)) oradius*sin(otest_pts(i) + (ostep/2))];
    for j = 1:2*N
        if j <= N
            %first object
            rnl = [-lambda + radius(1)*cos(test_pts(j) - (step/2)), radius(1)*sin(test_pts(j) - (step/2))];
            rnu = [-lambda + radius(1)*cos(test_pts(j) + (step/2)), radius(1)*sin(test_pts(j) + (step/2))];
            nhat= [-lambda + cos(test_pts(j)), sin(test_pts(j))];
        else
            %second object
            rnl = [lambda + radius(2)*cos(test_pts(j-N) - (step/2)), radius(2)*sin(test_pts(j-N) - (step/2))];
            rnu = [lambda + radius(2)*cos(test_pts(j-N) + (step/2)), radius(2)*sin(test_pts(j-N) + (step/2))];
            nhat= [lambda + cos(test_pts(j-N)), sin(test_pts(j-N))];
        end

        gk1  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k1);
        ggk1 = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k1);
        gk2  = @(m) green(rpnl,rpnu,0.5,rnl,rnu,m,k2);
        ggk2 = @(m) gradg(rpnl,rpnu,0.5,rnl,rnu,m,nhat,k2);

        if j<=N
            farfield(i) = farfield(i) - l(1) * (x(j) * gaussquad(gk1,n) - x(j+N) * gaussquad(ggk1,n));
        else
            farfield(i) = farfield(i) - l(2) * (x(j) * gaussquad(gk2,n) - x(j+N) * gaussquad(ggk2,n));
        end
    end
end

%% plotting the result

%plotting the results of VIE for verification
volume_two_disk

%SIE results
s = polarplot(onodes,-0.8*20*log10(2*pi*oradius*abs(farfield)),'blue');
set(s,'LineWidth',3);
legend({' VIE',' SIE'},'Location','northeast','Orientation','vertical')
ax = gca; 
ax.FontSize = 25; 
toc

%% defining functions

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
    [xi,wt] = lgwt(n,0,1);
    y = 0;
    for i = 1:n
        y = y + wt(i) * f(xi(i));
    end
end
