function ksp = Encode_low_memory(obj,p,traj)

Ns  = p.Ns;
Ne  = p.Ne;
Nc  = p.Nc;
Nx  = p.Nx_o;
Ny  = p.Ny_o;
Nsl = p.Nsl;
Nm  = p.Nm;

Nb0 = length(p.B0_e.vector);
Nt2 = length(p.T2_e.vector)/Nm;
Ntp = length(p.tp_e.vector);

%% Get B0 operator
if p.B0_e.flag == 1
    B0 = p.B0_e.val;
    
    ts = traj.ts;
    TE = p.TE;
    time = reshape(ts+TE,[Ns*Ne,1]);
    
    % load B0 map if file specified
    if ~isempty(p.B0_e.file)
        B0_map = [];
        load(p.B0_e.file) % Nx,Ny,Nsl
        
        % adapt number of slices to match desired Nsl
        if p.B0_e.Nsl ~= size(B0_map,3)
            if p.B0_e.Nsl == 1
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
    end

else
    B0 = 1;
end

%% loop over variables
for b = 1:Nb0 % loop over B0 offsets
    for t = 1:Nt2 % loop over B0 offsets
        for tp = 1:Ntp % loop over through-plane B0 gradients
            for sl = 1:Nsl % loop over slices
                SetEncodingValues();
                
                %% get operators except B0
                FO = GetFO(p.FO_e,p.res_o,p.FOV,traj.k);
                CS = GetCS(p.CS_e,p.fm,p.TE,traj.ts);
                T2 = GetT2(p.T2_e,p.TE,traj.ts);
                
                %% get B0 operator
                B0_s = B0_map(:,:,sl); % only use one slice for this iteration
                
                %% encode
                for d = p.d_range
                    x = obj(:,:,:,d);
                    x = x(:);
                        for s = 1:Ns
                            pos = s;
                            for e = 1:Ne-1
                                pos = [pos;s+(e*Ns)];
                            end
                        FO_tmp = FO(s,:,:,:);
                        if p.CS_e.flag == 1
                            CS_tmp = CS(s,:,:,:,:);
                        else
                            CS_tmp = 1;
                        end
                        if p.B0_e.flag == 1              
                            %% define field inhomogeneity operator
                            B0 = exp(-1i*2*pi*B0_s(:)*time(pos).');
                            try % B0 map is specified
                                B0_tmp = reshape(B0,[1,Ne,1,Nx*Ny,1]);
                            catch % only one value is specified
                                B0_tmp = reshape(B0,[1,Ne,1,1,1]);
                            end
                        else
                            B0_tmp = 1;
                        end
                        if p.T2_e.flag == 1
                            T2_tmp = T2(s,:,:,:,:);
                        else
                            T2_tmp = 1;
                        end
                        E = reshape(FO_tmp.*CS_tmp.*B0_tmp.*T2_tmp,[1*Ne*Nc,Nx*Ny*Nm]);
                        ksp(pos,d,b,t,tp,sl) = E*x;
                        end
                end
            end
        end
    end
end

% mean over all slices
ksp = sum(ksp,6)/Nsl;

function SetEncodingValues()
    p.tp_e.val = linspace(-p.tp_e.vector(tp)/2,p.tp_e.vector(tp)/2,Nsl);
    p.B0_e.val = p.B0_e.vector(b)+p.tp_e.val(sl);
    for m = 1:Nm
        p.T2_e.val(m) = p.T2_e.vector(Nm*(t-1)+m);
    end
end


end