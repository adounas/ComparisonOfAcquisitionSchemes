function ksp = Encode(obj,p,traj)

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

try % see if we run into memory issues

%% loop over variables
for b = 1:Nb0 % loop over B0 offsets
    for t = 1:Nt2 % loop over B0 offsets
        for tp = 1:Ntp % loop over through-plane B0 gradients
            for sl = 1:Nsl % loop over slices
                SetEncodingValues();
                
                %% get operators
                FO = GetFO(p.FO_e,p.res_o,p.FOV,traj.k);
                CS = GetCS(p.CS_e,p.fm,p.TE,traj.ts);
                try
                    B0 = GetB0(p.B0_e,p.res_o,p.FOV,p.TE,traj.ts,sl);
                    B0_memory_flag = 1;
                catch
                    B0_memory_flag = 0;
                end
                T2 = GetT2(p.T2_e,p.TE,traj.ts);
                
                %% encode
                for d = p.d_range
                    x = obj(:,:,:,d);
                    x = x(:);
                    try
                        E = FO.*CS.*B0.*T2;
                        E = reshape(E,[Ns*Ne*Nc,Nx*Ny*Nm]);
                        ksp(:,d,b,t,tp,sl) = E*x;
                    catch % if out of memory
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
                            if B0_memory_flag == 0
                                B0 = GetB0_low_memory(p.B0_e,p.res_o,p.FOV,p.TE,traj.ts,sl);
                            else
                                B0_tmp = B0(s,:,:,:);
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
end

catch
    ksp = Encode_low_memory(obj,p,traj);
    return
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