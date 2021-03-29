function [img,F_orig] = Decode_LCurve_iter(ksp,p,traj)

ksp = single(ksp);

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
Nit  = length(p.nIter.vector);

%% get operators
FO = single(GetFO(p.FO_d,p.res_i,p.FOV,traj.k));
CS = single(GetCS(p.CS_d,p.fm,p.TE,traj.ts));
B0 = single(GetB0(p.B0_d,p.res_i,p.FOV,p.TE,traj.ts));
T2 = 1;

%% decode
F = FO.*CS.*B0.*T2;
F = reshape(F,[Ns*Ne*Nc,Nx*Ny*Nm]);
F_orig = F;

img = zeros(Nx*Ny*Nm,Nd,Nb0,Nt2,Ntp,Nit);

for it = 1:Nit
    p.iter.val = p.iter.vector(it);
    F = F_orig;
    
    %% add L1 regularization if lambda is defined
%     if p.lambda.val ~= 0
%         F   = single([F  ; sqrt(p.lambda.val).*eye(size(F,2))]);
%     end
    
    for tp = 1:Ntp % loop over through-plane B0 gradients
        for b = 1:Nb0 % loop over all encoded B0 offsets
            for t = 1:Nt2 % loop over all encoded B0 offsets
                for d = p.d_range
                    
                    ksp_vec = single(ksp(:,d,b,t,tp));
%                     %% adapt ksp if lambda ~= 0
%                     if p.lambda.val ~= 0
%                         ksp_vec = [ksp_vec; zeros(size(F,2),1)];
%                     end
                    
                    switch p.solver
                        case 'CG'
                            rhs = F'*ksp_vec;
                            afun = @(x) F'*(F*x(:));
                            
                            [img_tmp,~,~,~,res] = pcg(afun, rhs, p.tol, p.iter.val, [], [], zeros(size(F,2),1));
                        case 'nlCG'
                            params.nIter = p.nIter;
                            img_tmp = nlCG(ksp_vec,F,[],params);
                    end
                    img(:,d,b,t,tp,it) = single(img_tmp);
                end
            end
        end
    end
end

img = reshape(img,[Nx,Ny,Nm,Nd,Nb0,Nt2,Ntp,Nit]);

end
