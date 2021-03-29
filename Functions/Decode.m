function img = Decode(ksp,p,traj)

Ns = p.Ns;
Ne = p.Ne;
Nc = p.Nc;
Nx = p.Nx_i;
Ny = p.Ny_i;
Nm = p.Nm;
Nd = p.Nd;

Nb0 = length(p.B0_e.vector);
Nt2 = length(p.T2_e.vector)/Nm;
Ntp = length(p.tp_e.vector);
Nl  = length(p.lambda.vector);

%% get operators
FO = GetFO(p.FO_d,p.res_i,p.FOV,traj.k);
CS = GetCS(p.CS_d,p.fm,p.TE,traj.ts);
B0 = GetB0(p.B0_d,p.res_i,p.FOV,p.TE,traj.ts);
T2 = 1;

%% decode
F = FO.*CS.*B0.*T2;
F = reshape(F,[Ns*Ne*Nc,Nx*Ny*Nm]);
F_orig = F;

img = zeros(Nx*Ny*Nm,Nd,Nb0,Nt2,Ntp,Nl);

for l = 1:Nl
    p.lambda.val = p.lambda.vector(l);
    F = F_orig;
    
    %% add L1 regularization if lambda is defined
    if p.lambda.val ~= 0
        F   = [F  ; sqrt(p.lambda.val).*eye(size(F,2))];
    end
    
    for tp = 1:Ntp % loop over through-plane B0 gradients
        for b = 1:Nb0 % loop over all encoded B0 offsets
            for t = 1:Nt2 % loop over all encoded B0 offsets
                for d = p.d_range
                    
                    ksp_vec = ksp(:,d,b,t,tp);
                    %% adapt ksp if lambda ~= 0
                    if p.lambda.val ~= 0
                        ksp_vec = [ksp_vec; zeros(size(F,2),1)];
                    end
                    
                    switch p.solver
                        case 'CG'
                            rhs = F'*ksp_vec;
                            afun = @(x) F'*(F*x(:));
                            
                            [img_tmp,~,~,~,res] = pcg(afun, rhs, p.tol, p.nIter, [], [], zeros(size(F,2),1));
                        case 'nlCG'
                            params.nIter = p.nIter;
                            img_tmp = nlCG(ksp_vec,F,[],params);
                    end
                    img(:,d,b,t,tp,l) = img_tmp;
                end
            end
        end
    end
end

img = reshape(img,[Ny,Nx,Nm,Nd,Nb0,Nt2,Ntp,Nl]);

end
