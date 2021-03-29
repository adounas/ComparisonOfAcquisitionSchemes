addpath(genpath('./Functions'));

p_list = [1,2];
for p_no = p_list

%% choose parameters
p = ChooseParameters(p_no);

%% load trajectory
traj = load(p.traj_file);
p.Ns = length(traj.ts);

%% load object
obj = LoadObject(p);

%% calculate encoding operator and encode
ksp = Encode(obj,p,traj);

%% add noise
ksp = AddNoise(ksp,p.std);

%% calculate decoding operator and decode
img = Decode(ksp,p,traj); % [Nx,Ny,Nm,Nd,Nb0,Nt2]

%% save results
save(p.save_file,'img','p','-v7.3');

end

clearvars -except img




