function obj = LoadObject(p)
% if object file is not defined, create a point object with Nsl slices
    if isempty(p.obj_file)
        obj = zeros(p.Nx_o,p.Ny_o,p.Nm,p.Nd,p.Nsl);
        for m = 1:p.Nm
            obj(round(size(obj,1)/2),round(size(obj,2)/2),m,m,:) = 1;
        end
    else
        obj = [];
        load(p.obj_file) % Nx,Ny,Nm,Nd,Nsl
        
        % only lactate for no chemical shift
        if p.CS_e.flag == 0
            obj = obj(:,:,1,:,:);
        end
        
        % adapt number of slices to match desired Nsl
        if p.Nsl ~= size(obj,5) 
            if p.Nsl == 1
                sl = round(size(obj,5)/2);
                obj = obj(:,:,:,:,sl);
            else
            sl = round(linspace(1,size(obj,5),p.Nsl))';
            sz = size(obj);
            tmp = zeros(sz(1),sz(2),sz(3),sz(4),length(sl));
            for i = 1:length(sl)
                tmp(:,:,:,:,i) = obj(:,:,:,:,sl(i));
            end
            obj = tmp;
            end
        end
        
    end
end