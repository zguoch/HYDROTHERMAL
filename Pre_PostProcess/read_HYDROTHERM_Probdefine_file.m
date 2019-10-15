function [data,varname] = read_HYDROTHERM_Probdefine_file(filename,nnod)

if nargin==0
    filename = 'Bench2/Out_pressure';
end

% Read node coordinates
fid = fopen(filename,'r');
data    = cell(1,4);

% Read porosity data
t          = skip_lines_until('--- Initial Porosity ---');
varname{1} = t;
% skip_lines(1);
t    = fgetl(fid);
i1 = strfind(t,'Uniform Value:')+14;
if (isempty(i1))
    data{1}    = read_data_blk(); % *NESTED FUNCTION*
else
    data{1}    = str2num(t(i1:end)); % *NESTED FUNCTION*
end



% Read permeability data
t          = skip_lines_until('--- X-Permeability ---',t);
varname{2} = t;
skip_lines(1);
data{2}    = read_data_blk(); % *NESTED FUNCTION*


% Read thermal conductivity data
t          = skip_lines_until('--- Thermal Conductivity ---',t);
varname{3} = t;
skip_lines(1);
data{3}    = read_data_blk(); % *NESTED FUNCTION*


% Read specific heat of rock data
t          = skip_lines_until('--- Specific Heat of Rock ---',t);
varname{4} = t;
% skip_lines(1);
t    = fgetl(fid);
i1 = strfind(t,'Uniform Value:')+14;
if (isempty(i1))
    data{4}    = read_data_blk(); % *NESTED FUNCTION*
else
    data{4}    = str2num(t(i1:end)); % *NESTED FUNCTION*
end



% Read rock density data
t          = skip_lines_until('--- Rock Density ---',t);
varname{5} = t;
% skip_lines(1);
t    = fgetl(fid);
i1 = strfind(t,'Uniform Value:')+14;
if (isempty(i1))
    data{5}    = read_data_blk(); % *NESTED FUNCTION*
else
    data{5}    = str2num(t(i1:end)); % *NESTED FUNCTION*
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

function data_blk = read_data_blk()
    data_blk = zeros(1,nnod);
    while ~feof(fid)
        t = fgetl(fid);
        if ~isempty(strfind(t,'---'))
            break
        end
        ind_x     = str2num(t);
        data_grid = [];
        while 1
            t     = fgetl(fid);
            if isempty(t)
                break
            end
            tmp   = str2num(t);
            ind_z = tmp(1);
            data_grid(ind_z,ind_x) = tmp(2:end);
        end
        skip_lines(2);
    end
    data_grid = data_grid(end:-1:1,:);
    data_grid = data_grid';
    data_blk  = data_grid(:);
end

end


