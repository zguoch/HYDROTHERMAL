function show_HYDROTHERM_1D(bench)

FigNo   = 100;
if nargin==0
    bench = '1_3cond';
end

Xrange  = [0 2000];
VarHTI  = load(['HTI_Sim_' bench]);
x_HTI   = VarHTI.x;
x_HTI   = x_HTI-x_HTI(1);
time    = VarHTI.time;

% Fluid properties
T_HTI   = VarHTI.Var.T;
P_HTI   = VarHTI.Var.P;
H_HTI   = VarHTI.Var.H;
Rho_HTI = VarHTI.Var.Rho2;
H_HTI   = 1e3 * H_HTI./Rho_HTI; % now in units of J/m^3
Mu_HTI  = VarHTI.Var.Mu1;
Vx      = VarHTI.Var.Vx;

% Matrix properties
Phi     = VarHTI.Var.Phi;
Perm    = VarHTI.Var.Perm;
K_m     = VarHTI.Var.K_m;
Rho_m   = VarHTI.Var.Rho_m;
Cp_m    = VarHTI.Var.Cp_m;

Trange  = [min(T_HTI(:))-10 max(T_HTI(:))+10];
Prange  = [min(P_HTI(:))-1e6 max(P_HTI(:))+1e6];

% 2) read min/max from HTI files; add 0.1*min/max; save ranges for benchmark

% PLOT
figure(FigNo);clf;
subplot(2,2,1);
    [ah,h1,h2] = plotyy(x_HTI,log10(Phi),x_HTI,log10(Perm));
    set(ah,'XLim',Xrange);
    set(h1,'Linestyle','--','LineWidth',2)
    set(h2,'Linestyle','-.','LineWidth',2);
    xlabel('Distance (m)');ylabel(ah(1),'log porosity');ylabel(ah(2),'log permeability (m^2)');
    
subplot(2,2,2);
    plot(x_HTI,Rho_m,'r-');
    set(gca,'XLim',Xrange);
    xlabel('Distance (m)');ylabel('Rock density (kg/m^3)');
        
subplot(2,2,3);
    plot(x_HTI,K_m,'r-');
    set(gca,'XLim',Xrange);
    xlabel('Distance (m)');ylabel('Thermal conductivity (W/m-K)');

subplot(2,2,4);
    plot(x_HTI,Cp_m,'r-');
    set(gca,'XLim',Xrange);
    xlabel('Distance (m)');ylabel('Rock density (J/kg-K)');

figure(FigNo+1);clf;
for ii=1:length(time)
    subplot(2,3,1);
        if ii>1; delete(h1); end
        h1 = plot(x_HTI,T_HTI(ii,:),'r-');
        title(sprintf(' Time = %0.4e s',time(ii)));
        if ii==1
            set(gca,'XLim',Xrange,'Ylim',Trange); hold all;
            xlabel('Distance (m)');ylabel('T (C)');
        end
    subplot(2,3,2);
        if ii>1; delete(h2); end
        h2=plot(x_HTI,P_HTI(ii,:),'r-');
        if ii==1
            set(gca,'XLim',Xrange,'Ylim',Prange); hold all;
            xlabel('Distance (m)');ylabel('P (Pa)');
        end
    subplot(2,3,4);
        if ii>1; delete(h3); end
        h3=plot(x_HTI,Rho_HTI(ii,:),'r-');
        if ii==1
            set(gca,'XLim',Xrange); hold all;
            xlabel('Distance (m)');ylabel('Densiy (kg/m3)');
        end
    subplot(2,3,5);
        if ii>1; delete(h4); end
        h4=plot(x_HTI,log10(Mu_HTI(ii,:)),'r-');
        if ii==1
            set(gca,'XLim',Xrange); hold all;
            xlabel('Distance (m)');ylabel('log Viscosity (Pa-s)');
        end
    subplot(2,3,3);
        if ii>1; delete(h5); end
        h5=plot(x_HTI,H_HTI(ii,:),'r-');
        title(sprintf(' Time = %0.4e d',time(ii)/(3600*24)));
        if ii==1
            set(gca,'XLim',Xrange); hold all;
            xlabel('Distance (m)');ylabel('Enthalpy (J/kg/K)');
        end
    subplot(2,3,6);
        if ii>1; delete(h6); end
        Htot_HTI = (1-Phi).*Rho_m.*Cp_m.*T_HTI(ii,:) + ...
           Phi.*Rho_HTI(ii,:).*H_HTI(ii,:);
        h6=plot(x_HTI,Htot_HTI,'r-');
        if ii==1
            set(gca,'XLim',Xrange); hold all;
            xlabel('Distance (m)');ylabel('Total Enthalpy');
        end
    drawnow
    pause(0.1);
end

end
