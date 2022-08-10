function with_cond_all_zloop_sinefn(z, Lx, tmax, dx, dt, startind, period, w0, Df, Ds, dd)

cd('test_for_github') %or where the results will be put
dirpath  = ['Figures_z_', num2str(z), '_Lx_' num2str(Lx) '_w0_', num2str(w0), '_sine_', num2str(period),'_d_',num2str(dd)];
mkdir(dirpath);

epnorm   = 1.5e-5; % break out of time loop when (Tf) is smaller than this (negative temps) - empirical
dTnorm   = 1;
nstreams = 10; % level of decimation when plotting
nout     = 50; % number of files that are written for each case

x    = [0:dx:Lx];
n    = length(x);
ninc = floor(n/nstreams);

vv   = 1./z;

% initialize
%Df = 0; only debugging
%Ds = 0; only debugging

T_solid = zeros(1, n); 
T_fluid = zeros(1, n);
T_diff  = zeros(1, n);

T_f_try = T_fluid;
T_diff  = zeros(1, n);
T_solid_steps = [];
T_fluid_steps = [];
t_steps       = [];
legendtext    = [];

dT_fluid_dx   = gradient(T_fluid, dx); % starting gradient in fluid
delta_T       = T_fluid - T_solid;     % starting contrast
diff_f        = Df*4*del2(T_fluid,dx); % KK*(starting 2nd derivative in fluid)
diff_s        = Ds*4*del2(T_solid,dx);

%Value of initial temperature prior to perturbation
T_in   = 0;
% Value of temp of fluid at left edge
T_0    = 0; %initial inlet temp

h1 = figure(1); grid on; box on; hold on
%h1 = subplot(211); grid on; box on; hold on
%h2 = subplot(212); grid on; box on; hold on
solid_plot = plot(x, T_solid);
fluid_plot = plot(x, T_fluid);
set(solid_plot, 'linewidth',[1]);
set(fluid_plot, 'linewidth',[1]);
ylim([-0.1,1]);
ylabel('Temperature/{\Delta}T'); xlabel('x');
set(gca,'fontname','helvetica','fontsize', [14])

diff_plot     = plot(x, T_diff);
[mdiff, indm] = max(T_diff);
mx(1,:)       = [0, x(indm)];
lm            = [mx 0; mx T_0];
max_plot      = plot(lm(:,1), lm(:,2));
set(max_plot,'linestyle','--','color','k')
legend('{T''}_s', '{T''}_f','{T''}_f - {T''}_s');

k    = 1;
step = 0;
maxiter = length([dt:dt:tmax]);
outiter = floor(maxiter/nout);
iter    = 1;
endflag = 0;
brflag  = 0;

for iter = 1:maxiter
    step  = step + 1;
    t     = dt*iter;
    
    % take an intermediate step to evolve the fluid temp 
    dtsub = dt*0.5;

    dT_fluid_dt = - delta_T - dT_fluid_dx + diff_f;
    T_f_try = T_fluid + vv * dT_fluid_dt * dtsub;
    % evolve solid temp at int step
    T_solid_try = T_solid + (delta_T + diff_s) * dtsub;
    % find fluid gradient at substep
    dT_fluid_dx_try = gradient(T_f_try, dx);
    delta_T         = T_f_try - T_solid_try;
    
    % use values after intermediate step in order to take next step
    dT_fluid_dt = - delta_T - dT_fluid_dx_try + diff_f;
    
    % take full step
    T_fluid = T_fluid + vv * dT_fluid_dt * dt;
    T_solid = T_solid + (delta_T + diff_s) * dt;
    
    %T_fluid(in) = 0.15*T_f_try(in)+(1-0.15)*T_fluid(in);
    delta_T     = T_fluid - T_solid;
    dT_fluid_dx = gradient(T_fluid, dx);


    diff_f        = Df*4*del2(T_fluid,dx); % KK*(2nd derivative in fluid)
    diff_s        = Ds*4*del2(T_solid,dx);
   
    %T_fluid(1)  = 0.5*(tanh((t-0.5*period)/w0) - tanh((t-1.5*period)/w0)); %0.5*(1+sin(2*pi*t/period + pi/2));  
    T_fluid(1)  = 0.5*(1+sin(2*pi*t/period + pi/2));
    % boundary condition on Tf only at left
    if iter == 1 
        T_fluid(2) = 0.5*T_fluid(1);
    end
    
    T_diff = T_fluid - T_solid;
    dTnorm = norm(T_fluid); 
    
    if (mod(iter, outiter) == 0)
        set(solid_plot, 'YDATA', T_solid);
        set(fluid_plot, 'YDATA', T_fluid);
        set(diff_plot,  'YDATA', T_diff);
        [mdiff, indm] = max(T_diff);
        mx        = [mx; t x(indm)];
        lm        = [x(indm) 0; x(indm) T_0];
        set(max_plot, 'XDATA', lm(:,1), 'YDATA', lm(:,2));
        drawnow
        tl = ['t = ' num2str(t)];
        title(tl);
        fn = [dirpath '/Temp_x_' num2str(t) '.png'];
        %print(fn, '-dpng')
	    fn = [dirpath '/Temp_x_' num2str(t) '.fig'];
	    savefig(fn)
        
        % save a growing array of all solutions through time
        Tf(k,:) = T_fluid;
        Ts(k,:) = T_solid;
        k       = k + 1;
    end
    
    % nstreams sets the level of decimation
    if mod(step, floor(1 / dt / nstreams)) == 0
        t;
        t_steps = [t_steps; t];
        % make a decimated sampling of Temperatures in the domain
        T_solid_steps = [T_solid_steps; T_solid(1:ninc:end)];
        T_fluid_steps = [T_fluid_steps; T_fluid(1:ninc:end)];
    end
    
   if sum(T_fluid<-epnorm) | sum(T_solid<-epnorm) > 0
       brflag = 1
       break
   elseif max(T_fluid>1.5) | max(T_solid>1.5) > 0
       brflag = 2
       break
   end
    
end

if brflag == 1
        %for debugging:
        %set(solid_plot, 'YDATA', T_solid);
        %set(fluid_plot, 'YDATA', T_fluid);
        %set(diff_plot,  'YDATA', T_diff);
    disp(['negative Temperature - unstable!, Df = ', num2str(Df), ' Ds = ', num2str(Ds), ' d = ', num2str(dd), ' z = ', num2str(z)])
    %iter;
    %[T_fluid' T_solid'];
    %pause;
elseif brflag == 2
    disp(['Temperature > 1 - unstable!, Df = ', num2str(Df), ' Ds = ', num2str(Ds), ' d = ', num2str(dd), ' z = ', num2str(z)])
    %iter;
    %[T_fluid' T_solid'];
else
    x_points = x(1:ninc:end); % make a decimated sampling of domain
    time_of_contact = x_points / vv;
    y = time_of_contact;
    legendtext = strcat('x''=',string(num2cell(x_points)));
    zz = (t_steps - time_of_contact); 
    
    h2=figure(); hold on; grid on; box on; title('Solid');
    xlabel('t'' = (t - x*z)'); ylabel('Ts');
    
    for ind = 1:length(y)
        plot(zz(:,ind), T_solid_steps(:,ind))
    end
    xlim([0,zz(end,1)]);
    axP = get(gca,'Position');
    legend(legendtext,'Location','NorthEastOutside')
    axP = get(gca,'Position');
    
    h3=figure(); hold on; grid on; box on; title('Fluid');
    xlabel('t'' = (t - x*z)'); ylabel('Tf');
    for ind = 1:length(y)
        plot(zz(:,ind), T_fluid_steps(:,ind))
    end
    xlim([0,zz(end,1)]);
    legend(legendtext,'Location','NorthEastOutside')
    
    h4=figure();
    hold on; grid on; box on; title('Location of Max(T_f - T_s) vs t');
    plot(mx(:,1),mx(:,2),'k.');
    %plot(mx(:,1),mx(:,2),'-');
    p = polyfit(mx(startind:end,1), mx(startind:end,2), 1);
    f = polyval(p, mx(:,1));
    plot(mx(:,1),f,'k--');
    xlabel('t'), ylabel('Location of max(T_f - T_s)')
    title(['Slope = ' num2str(p(1)) ])
    
    figure(2); set(gca,'fontname','helvetica','fontsize', [14])
    f1 = [dirpath '/Tsolid.png'];
    %print(f1, '-dpng')
    f1 = [dirpath '/Tsolid.fig'];
    savefig(f1);
    figure(3); set(gca,'fontname','helvetica','fontsize', [14])
    f2 = [dirpath '/Tfluid.png'];
    %print(f2, '-dpng')
    f2 = [dirpath '/Tfluid.fig'];
    savefig(f2);
    figure(4); set(gca,'fontname','helvetica','fontsize', [14])
    f3 = [dirpath '/Maxt.png'];
    %print(f3, '-dpng')
    f3 = [dirpath '/Maxt.fig'];
    savefig(f3);
    f4 = [dirpath '/vars.mat'];
    save(f4); % write all vars to a .mat file for later use in analysis

end

cd ../


end

