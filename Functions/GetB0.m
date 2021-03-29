function B0 = GetB0(opts,res,FOV,TE,ts,slice)

if opts.flag == 1
    B0 = opts.val;
    Ns = length(ts);
    Ne = length(TE);
    Nx = round(FOV(1)/res(1));
    Ny = round(FOV(2)/res(2));
    
    t = reshape(ts+TE,[Ns*Ne,1]);
    
    % load B0 map if file specified
    if ~isempty(opts.file)
        B0_map = [];
        load(opts.file) % Nx,Ny,Nsl
        
        % adapt number of slices to match desired Nsl
        if opts.Nsl ~= size(B0_map,3)
            if opts.Nsl == 1
                sl = round(size(B0_map,3)/2);
                B0_map = B0_map(:,:,sl);
            else
                sl = round(linspace(1,size(B0_map,3),p.Nsl))';
                sz = size(B0_map);
                tmp = zeros(sz(1),sz(2),length(sl));
                for i = 1:length(sl)
                    tmp(:,:,:,:,i) = B0_map(:,:,sl(i));
                end
                B0_map = tmp;
            end
        end
        B0 = B0_map(:,:,slice); % only use one slice for this iteration
    end

    %% define field inhomogeneity operator
    B0 = exp(-1i*2*pi*B0(:)*t.');
    try % B0 map is specified
        B0 = permute(B0,[2,1]);
        B0 = reshape(B0,[Ns,Ne,1,Nx*Ny,1]);
    catch % only one value is specified
        B0 = reshape(B0,[Ns,Ne,1,1,1]);
    end
else
    B0 = 1;
end

end




