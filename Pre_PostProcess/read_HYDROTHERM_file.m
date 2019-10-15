function [data,istep,time,varnames,x,y,z,nx,ny,nz] = read_HYDROTHERM_file(filename)

if nargin==0
    filename = 'Bench2/Out_pressure';
end

% Read node coordinates
fid = fopen(filename,'r');
skip_lines_until('X-Direction'); skip_lines(1);
x   = []; nblk = 0;
while ~feof(fid)
    t       = fgetl(fid);
    if isempty(t)
        break
    end
    inod    = str2num(t);
    t       = fgetl(fid);
    x(inod) = str2num(t);
    skip_lines(3);
    nblk = nblk + 1;
end
skip_lines_until('Y-Direction'); skip_lines(1);
nx  = length(x);

y   = [];
while ~feof(fid)
    t       = fgetl(fid);
    if isempty(t)
        break
    end
    inod    = str2num(t);
    t       = fgetl(fid);
    y(inod) = str2num(t);
    skip_lines(3);
end
ny  = length(y);

skip_lines_until('Z-Direction'); skip_lines(1);
z   = [];
while ~feof(fid)
    t       = fgetl(fid);
    if isempty(t) || ~isempty(strfind(t,'Time'))
        break
    end
    inod    = str2num(t);
    t       = fgetl(fid);
    z(inod) = str2num(t);
    skip_lines(3);
end
nz  = length(z);

if length(y)==1
    if length(z)>1
        [x2,z2] = meshgrid(x,z(end:-1:1));
        x2 = x2';
        x2 = x2(:)';
        z2 = z2';
        z2 = z2(:)';
        x  = x2;
        z  = z2;
    end
end

nnod = length(x);

% Read time step data
time = []; iout = 0; ivar = 0;
data = cell(1);

t = skip_lines_until('Time Step No.',t);
while ~feof(fid)
    iout        = iout + 1;
    
    % Read time step number
    i1 = strfind(t,'Time Step No.')+13;
    i2 = strfind(t,'; Simulation time:')-1;
    istep(iout) = str2num(t(i1:i2));
    % Read current simulation time
    i1 = i2+19;
    i2 = strfind(t,'(s)')-1;
    if(isempty(i2))  %if (s) doesn't exist, try (yr)
        i2 = strfind(t,'(yr)')-1;
    end
    time(iout)  = str2num(t(i1:i2));
    
    skip_lines(1);
    
    t    = fgetl(fid);
    while ~feof(fid)
        ivar               = ivar + 1;
        data{ivar}(:,iout) = zeros(nnod,1);
        varnames{ivar}      = t;
        
        t = fgetl(fid);
        if strfind(t,'Maximum value')
            skip_lines(1);
        end
        
        data_grid = zeros(nz,nx);
        for iblk=1:nblk
            t = fgetl(fid);
            if strcmp(t,'++ Row  1 ++')
                if iout==1
                    data{ivar} = fscanf(fid,'%f');
                end
                break
            elseif isempty(t)
                break
            else
                ind_x = str2num(t);
            end
            
            for iz = 1:nz
                t     = fgetl(fid);
                tmp   = str2num(t);
                ind_z = tmp(1);
                data_grid(ind_z,ind_x) = tmp(2:end);
            end
            
            skip_lines(2);
            t       = fgetl(fid);
        end
        
        data_grid  = data_grid(nz:-1:1,:);
        data_grid  = data_grid';
        data{ivar}(:,iout) = data_grid(:);
        
        while ~feof(fid)
            t = fgetl(fid);
            if strfind(t,'Time Step No.')
                new_time_step=1;
                break
            elseif strfind(t,'---')
                new_time_step=0;
                break
            end
        end
        
        if new_time_step
            ivar = 0;
            break 
        end
    end
    % t=skip_lines_until(t,'Time Step No.',t);
end
fclose(fid);


% -------------------------- NESTED FUNCTIONS 


function skip_lines(n)
    for i=1:n
        fgetl(fid);
    end
end

function t=skip_lines_until(string,t)
    if nargin>1 && ~isempty(t) && ~isempty(strfind(t,string))
        return
    end
    while ~feof(fid)
        t = fgetl(fid);
        if ~isempty(strfind(t,string))
            break
        end
    end
end

end


