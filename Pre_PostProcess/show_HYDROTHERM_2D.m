function show_HYDROTHERM_2D()

FigNo    = 100;
folder = '../examples/Ex.09';
file     = 'HTI_data.mat';
var2plot = {'V_l' 'V_v' 'T' 'P' 'H' 'Mu_l' 'Mu_v' 'Rho_l' 'Rho_v'};
make_png = 1;

DATA     = load([folder '/' file]);
nx       = DATA.nx;
nz       = DATA.nz;
xg       = reshape(DATA.x,nx,nz);
zg       = reshape(DATA.z,nx,nz);
times    = DATA.time;
VAR      = DATA.VAR;

sec_2_yr = 1/(365.25*24*3600);

for ivar=1:length(var2plot)
    varname   = var2plot{ivar};
    switch varname
        case 'V_l'
            var = sqrt(VAR.Vx_l.^2+VAR.Vz_l.^2);
        case 'V_v'
            var = sqrt(VAR.Vx_v.^2+VAR.Vz_v.^2);
        otherwise
            var = VAR.(varname)
            vec = [];
    end
    
    clims = [min(min(var)) max(max(var))];
    for ii=1:length(times)
        cg  = reshape(var(ii,:),nx,nz);
        figure(FigNo);clf
        contourf(xg,zg,cg,100,'Linecolor','none'); hold on
        axis equal tight
        colorbar
        title(sprintf('%s at t=%.2f yr',varname,times(ii)*sec_2_yr));
        if min(clims)<max(clims)
            set(gca,'CLim',clims);
        end
        
        switch varname
            case 'V_l'
                quiver(xg(:)',zg(:)',VAR.Vx_l(ii,:),VAR.Vz_l(ii,:),1,'k');
            case 'V_v'
                quiver(xg(:)',zg(:)',VAR.Vx_v(ii,:),VAR.Vz_v(ii,:),1,'k');
        end

        if make_png
            filename = ['HTI_' varname '_' num2str_d(ii,4)];
            print('-dpng','-r150',[folder '/' filename]);
        end
    end
end

end