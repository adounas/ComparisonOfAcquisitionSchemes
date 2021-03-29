function FO = GetFO(opts,res,FOV,k)

if opts.flag == 1
    Ns = length(k);
    
    %% define grid
    Nx = round(FOV(1)/res(1));
    Ny = round(FOV(2)/res(2));
    
    dx = FOV(1)/Nx;
    dy = FOV(2)/Ny;
    [xq,yq] = meshgrid((-floor(Nx/2):ceil(Nx/2)-1)*dx,...
                       (-floor(Ny/2):ceil(Ny/2)-1)*dy);
    grid = [reshape(xq,[Nx*Ny,1]),reshape(yq,[Nx*Ny,1])];
    %% define Fourier operator
    FO = 1/(Nx*Ny)*exp(-1i*k*grid.');
    FO = reshape(FO,[Ns,1,1,Nx*Ny,1]);
else
    FO = 1;
end

end