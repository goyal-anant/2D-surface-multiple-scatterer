function volume_disk
    %%% defining constants
    % epsilon0 = 8.853e-12;
    epsilonr = 4;           %relative permitivity of the dielectric shell
    epsilon  = epsilonr;     
    theta0   = pi/2;        %angle between z-axis and propagation axis
    phi0     = pi;           %angle between x-axis and projection of k on xy plane
    lambda   = 1;           %wavelength 
    E0       = 1;           %electric field amplitude 
    k        = 2*pi/lambda; %wave vector
    outer    = 0.30;        %outer radius
    inner    = 0.25;        %inner radius
    
    %%%building the scatterer
    %function [X1, Y1] = dielectric_shell(outer,inner,lambda,N,plot_flag)
    N         = 50; %no. of pixels in the space
    plot_flag = 0;  %put 1 to plot
    % [X1,Y1,r,X,Y]   = dielectric_shell(outer,inner,lambda,N,plot_flag);
    % [X1, Y1, r, X, Y] = dielectric_half_shell(outer,inner,lambda,N,plot_flag);
    [X1, Y1] = structure(N,lambda,outer,plot_flag);
    
    %%%formulating and solving the problem
    %a = %lambda/50;  %side of the square patch
    r = (2*lambda/N)/sqrt(pi); %radius of equivalent circle with same cross section
    
    %forming matrix C in the equation [C][E] = Ei
    for n = 1:length(X1)
        for m = 1:length(X1)
            if n == m
                C(n,m) = 1 + ((epsilon - 1) * (1i/2) * (pi * k * r * besselh(1,2,(k * r)) - 2i));
            else
                rho_mn = sqrt((X1(n) - X1(m))^2 + (Y1(n) - Y1(m))^2);
                C(n,m) = (1i * pi * k * r / 2) * (epsilon - 1) * (besselj(1, (k * r))) * (besselh(0,2,(k * rho_mn)));
            end
        end
    end
    
    %for Ei vector
    alpha = (X1 * sin(theta0) * cos(phi0)) + (Y1 * sin(theta0) * sin(phi0));
    Ei = E0 * exp(-1i * k * alpha);     %incident electric field 
    
    %solving matrix equation
    E = C\Ei';      %total electric field on the dielectric volume
    
    %%%getting the plot of E in terms of phi
    % temp = 1;
    % for i = 1:length(X1)
    %     if X1(i) >= 0 && Y1(i) >= 0
    %         phi(temp) =  0 + abs(atan(Y1(i)/X1(i)));
    %         temp = temp + 1;
    %     elseif X1(i) <= 0 && Y1(i) >= 0
    %         phi(temp) =  pi - abs(atan(Y1(i)/X1(i)));
    %         temp = temp + 1;
    %     elseif X1(i) <= 0 && Y1(i) <= 0
    %         phi(temp) =  pi + abs(atan(Y1(i)/X1(i)));
    %         temp = temp + 1;
    %     elseif X1(i) >= 0 && Y1(i) <= 0
    %         phi(temp) =  2*pi - abs(atan(Y1(i)/X1(i)));
    %         temp = temp + 1;
    %     end
    % end
    % 
    % phi = sort(phi)*180/pi;
    % It does not matter where you start the phi from, whereever you start it
    % will keep on repeating, so generate a random phi from 0 to 200 to match
    % with richmond's plot and plot abs(E) wrt the phi generated
    
    % phi_E = linspace(0,200,length(E));
    % AxesH = axes('YTick',0:0.1:2, 'NextPlot', 'add');
    % plot(phi_E,flip(abs(E)),linewidth = 3); 
    % hold on; grid on; xlabel('\phi(degrees)'); ylabel('|E|');
    % set(gca,'fontsize',20)
    % 
    % import_rr = load('extrd_fig3_rich.mat');
    % rr = import_rr.rr;
    % plot(rr(:,1),rr(:,2),linewidth = 3);
    
    % error_vec = abs(E) - rr(:,2);
    % error = norm(error_vec);
    
    % now the problem is the wiggles in the plot, the envelope is matching with
    % the results of figure 3 in richmond's paper. Why are the wiggles coming?
    % What is a possible solution.
    
    
    %the problem I am facing here is that E contains the total electric field
    %in cell 1,2,3,... and I do not have the coordinate mapping of cells and
    %(X1,Y1). So, when I am plotting E wrt anything on x-axis, it is just
    %showing the values of E stored in cell 1,2,3,... irrespective of the
    %x-axis
    
    
    %%%plotting the echowidth
    %distant scattered field
    phi = linspace(0,2*pi,length(E));   %for circular variation around the shell
    
    %%%for distant scattered field
    rho0 = 500*lambda;          %distance where we want to check the field
    %temporary variable 'temp' to store the amplitude of the field
    temp = -1i * (pi * k / 2) * sqrt(2i / (pi * k * rho0)) ...
        * exp(-1i * k * rho0) * (epsilon-1) * r * besselj(1,(k * r)) ;
    temp1 = cos(phi.')*X1 + sin(phi.')*Y1;  % to hold angular variation
    Es = temp * exp(1i * k * temp1) * E;  % distant scattered field
    
    
    %%%Wphi
    %temporary variable 'temp' to store the amplitude of the field
    % temp = -1i * (pi * k / 2) * sqrt(2i / (pi * k * rho0)) ...
    %     * exp(-1i * k * rho0) * (epsilon-1) * r * besselj(1,(k * r)) ;
    temp1 = cos(phi.')*X1 + sin(phi.')*Y1;  % to hold angular variation
    % Es = temp * exp(1i * k * temp1) * E;  % distant scattered field
    
    t1 = (epsilon-1) * r * besselj(1,(k*r));
    t2 = t1*exp(1i * k * temp1) * E;
    t2 = abs(t2).^2;
    t3 = 1;%abs(Ei).^2;
    
    Wphi = (k*pi^2) * t2 / t3';
    
    %%%plotting rcs(dB)
    % AxesH = axes('XTick',0:10:200, 'NextPlot', 'add');
    % plot(phi*180/pi,Wphi,'linewidth',3);
    v = polarplot(phi,20*log10(2*pi*rho0*abs(Es)),'--r');
    set(v,'LineWidth',3)
    hold on;
    % hold on; grid on; set(gca,'fontsize',20)
    % xlabel('\phi(degrees)');
    % ylabel('Echo width/ Wavelength');
end

