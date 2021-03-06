function p = ChooseParameters(p_no)
%% choose parameters for subsequent simulation
% definitions:
%   FOV     [m,m]: field-of-view in x and y
%   res_o   [m,m]: object resolution in x and y
%   res_i   [m,m]: image resolution in x and y
%   fm    Nm*[Hz]: chemical shift

% initialize parameter struct
Initialize();

%% case summary:
% F: Fourier operator
% C: chemical shift
% B: field inhomogeneities
% T: T2* relaxation
% t: through-plane B0 gradients
% i: XCAT object (if not then point object)
% b: loaded B0 map
% n: noise
% r: regularization
% p: partial voluming in x/y
% l: large FOV

% x=1: optimal trajectory
% x=2: EPI
% x=3: SPI
%        F C B T t i b n r p l - point object and operators
% x01:  |x| | | | | | | | | | |
% x02:  |x|x| | | | | | | | | |
% x03:  |x| |x| | | | | | | | |
% x04:  |x| | |x| | | | | | | |
% x05:  |x|x|x| | | | | | | | |
% x06:  |x|x| |x| | | | | | | |
% x07:  |x| |x|x| | | | | | | |
% x08:  |x|x|x|x| | | | | | | |
% x09:  |x| |x| |x| | | | | | |
% x10:  |x|x|x|x|x| | | | | | |

%        F C B T t i b n r p l - XCAT and B0 map
% x11:  |x| | | | |x| | | | | |
% x12:  |x|x| | | |x| | | | | |
% x13:  |x| |x| | |x| | | | | |
% x14:  |x| |x| |x|x| | | | | |
% x15:  |x| |x| | |x|x| | | | |
% x16:  |x| |x| |x|x|x| | | | |

%        F C B T t i b n r p l - noise, regularization, partial voluming
% x21:  |x| | | | |x| |x| | | |
% x22:  |x| | | | |x| |x|x| | | % cf. special case 100x
% x23:  |x| | | | | | | | |x| |
% x24:  |x|x| | | | | | | |x| |
% x25:  |x| | | | |x| | | |x| |
% x26:  |x|x| | | |x| | | |x| |
% x27:  |x|x|x|x|x|x| | | |x| |

%        F C B T t i b n r p l - large FOV - long calculations
% x31:  |x| | | | | | | | | |x|
% x32:  |x| | | | | | | | |x|x|
% x33:  |x|x|x|x|x|x| | | | |x|
% x34:  |x|x|x|x|x|x| | | |x|x|
% x35:  |x|x|x|x| |x|x| | | |x|
% x36:  |x|x|x|x| |x|x| | |x|x|
% x37:  |x|x|x|x| |x|x|x|x| |x|
% x38:  |x|x|x|x| |x|x|x|x|x|x|

% special cases:
% 1001: Calculate L curve for small FOV - use code "CoAS_Sim_GetLCurve.m" -
% optimal traj
% 2001: Calculate L curve for large FOV with pv (OT) - use code "CoAS_Sim_GetLCurve.m"
% 2002: Calculate L curve for large FOV with pv (SPI)- use code "CoAS_Sim_GetLCurve.m"
% 2003: Calculate L curve for small FOV with pv (EPI)- use code "CoAS_Sim_GetLCurve.m"

% for paper:
% x0'001: B0 sweep, large FOV, pv, PSF, no reg
% x0'002: T2 sweep, large FOV, pv, PSF, no reg
% x0'003: tpB0 sweep, large FOV, pv, PSF, no reg
% x0'004: B0 sweep, large FOV, pv, PSF, reg
% x0'005: T2 sweep, large FOV, pv, PSF, reg
% x0'006: tpB0 sweep, large FOV, pv, PSF, reg

% x0'011: XCAT, no pv, no B0, no T2*, no noise, no reg
% x0'012: XCAT, pv, no B0, no T2*, no noise, no reg
% x0'013: XCAT, pv, B0, T2*, no noise, no reg
% x0'014: XCAT, pv, B0, T2*, noise, reg

% 100'001-100'002 noise test for code CoAS_Sim_noise_test

% cases for Valery
% 1: EPI as for in vivo - high/low B0, T2*
% 3: SPI as for in vivo - high/low B0, T2*

%% case switch
switch p_no
       case 1
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.025];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,20];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed

        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
       case 2
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.025];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,20];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed

        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
         case 2005 % calc large FOV L curve with pv - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/EPI_traj_201210.mat';

        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 2
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.lambda.vector = [0,0.0001];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
        
    case 101 % only FO, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.save_file = ['./Results/',num2str(p_no),'.mat'];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 102 % FO,CS, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 103 % FO,B0, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 104 % FO,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 105 % FO,CS,B0, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 106 % FO,CS,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 107 % FO,B0,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 108 % FO, CS, B0, and T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 109 % FO,B0,tpB0, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 110 % FO, CS, B0, and T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 111 % FO, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 112 % FO,CS XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 113 % FO,B0 map, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 114 % FO,B0 map, XCAT, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 115 % FO,B0 map, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20_nopv.mat';
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 116 % FO,B0 map, XCAT, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20_nopv.mat';
        p.Nsl = 9; % tp does not need to be specified for loaded B0 map
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 121 % FO, XCAT, no partial voluming, small FOV, noise
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';

        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 122 % FO, XCAT, no partial voluming, small FOV, noise, reg
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';

        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        p.lambda.vector = [0,0.01];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 123 % FO, XCAT, no partial voluming, small FOV, noise
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 124 % FO, XCAT, no partial voluming, small FOV, noise
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 125 % FO, XCAT, no partial voluming, small FOV, noise
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';

        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 126 
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 127 
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 131 % only FO, no partial voluming
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        p.save_file = ['./Results/',num2str(p_no),'.mat'];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 132 % only FO
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        p.save_file = ['./Results/',num2str(p_no),'.mat'];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 133 
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 134
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 135
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
             
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40_nopv.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 136
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 137
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
             
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40_nopv.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        p.lambda.vector = [0,0.01];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 138
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.01];
        
        p.traj_file = './Data/optimal_trajectory.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        p.lambda.vector = [0,0.01];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;

        
     case 301 % only FO, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 302 % FO,CS, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 303 % FO,B0, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 304 % FO,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 305 % FO,CS,B0, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 306 % FO,CS,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 307 % FO,B0,T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 308 % FO, CS, B0, and T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 309 % FO,B0,tpB0, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 310 % FO, CS, B0, and T2, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 50;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [Inf,0.001];
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 311 % FO, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 312 % FO,CS XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 313 % FO,B0 map, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 314 % FO,B0 map, XCAT, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.vector = [0,100];
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.tp_e.flag = 1;
        p.tp_e.vector = [0,50];
        p.Nsl = 9;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 315 % FO,B0 map, XCAT, no partial voluming, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20_nopv.mat';
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 316 % FO,B0 map, XCAT, no partial voluming in x/y, small FOV
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20_nopv.mat';
        p.Nsl = 9; % tp does not need to be specified for loaded B0 map
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
      
     case 1001 % calc small FOV L curve - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory_halfFOV.mat';

        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 2001 % calc large FOV L curve with pv - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/optimal_trajectory.mat';

        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 2002 % calc large FOV L curve with pv - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210.mat';

        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 2003 % calc large FOV L curve with pv - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/EPI_traj_201210.mat';

        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
     case 1003 % calculate L curve - use with code "CoAS_Sim_GetLCurve.m"
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210_halfFOV.mat';

        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 17; % range to be reconstructed
        
        p.std = 0.002;
        l = 1:0.2:5;
        lambda = 10.^-l;
        p.lambda.vector = [lambda,0];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
        
    case 20001
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = -50:5:50;
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 20002
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [5:5:100,Inf]/1000;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 20003
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.tp_e.flag = 1;
        p.tp_e.vector = 0:5:100;
        p.Nsl = 9;
        p.Nsl = 99;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 20004
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = -50:5:50;
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.lambda.vector = 0.001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 20005
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [5:5:100,Inf]/1000;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.lambda.vector = 0.001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 20006
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.tp_e.flag = 1;
        p.tp_e.vector = 0:5:100;
        p.Nsl = 9;
        p.Nsl = 99;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.lambda.vector = 0.001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30001
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = -50:5:50;
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30002
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [5:5:100,Inf]/1000;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30003
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.tp_e.flag = 1;
        p.tp_e.vector = 0:5:100;
        p.Nsl = 9;
        p.Nsl = 99;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30004
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.B0_e.vector = -50:5:50;
        p.B0_d.flag = 0;
        p.B0_d.val = 0;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.lambda.vector = 0.0001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30005
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.T2_e.flag = 1;
        p.T2_e.vector = [5:5:100,Inf]/1000;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.lambda.vector = 0.0001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
    case 30006
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.B0_e.flag = 1;
        p.tp_e.flag = 1;
        p.tp_e.vector = 0:5:100;
        p.Nsl = 9;
        p.Nsl = 99;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.lambda.vector = 0.0001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 20011
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_20_20_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed

        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 20012
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 20013
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = 0.025;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 20014
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        %p.TE0 = 0.0015855;
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_20_20.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = 0.025;
        
        p.traj_file = './Data/EPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_20_20.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.std = 0.002;
        p.lambda.vector = 0.0001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 30011
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.005, 0.005];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_40_40_nopv.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 30012
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 30013
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = 0.025;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 30014
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0015855;
        p.dTE = 0.0011500;
        p.Ne  = 6;
           
        p.B0_e.flag = 1;
        p.B0_e.file = './Data/B0_40_40.mat';
        
        p.T2_e.flag = 1;
        p.T2_e.vector = 0.025;
        
        p.traj_file = './Data/SPI_traj_201210.mat';
        
        p.obj_file  = './Data/XCAT_40_40.mat';
        p.Nd = 40;      % number of dynamics in obj_file
        p.d_range = 1:40; % range to be reconstructed
        
        p.Nsl = 9;
        
        p.std = 0.002;
        p.lambda.vector = 0.0001;
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 100001 % noise test EPI for CoAS_Sim_noise_test
        p.FOV   = [0.200, 0.200];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/EPI_traj_201210.mat';

        p.std = 0.002;
        p.lambda.vector = [0,0.001];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
   case 100002 % noise test SPI for CoAS_Sim_noise_test
        p.FOV   = [0.400, 0.400];
        p.res_o = [0.001, 0.001];
        p.res_i = [0.005, 0.005];
        
        p.CS_e.flag = 1;
        p.CS_d.flag = 1;
        p.met = {'lac','pyH','ala','pyr','bic'};
        p.fm  = [-388.7853, -273.1137, -167.0813, 0, 330.9495];
        p.TE0 = 0.0096735;
        p.dTE = 0.0011500;
        p.Ne  = 6;
        
        p.traj_file = './Data/SPI_traj_201210.mat';

        p.std = 0.002;
        p.lambda.vector = [0,0.0001];
        
        p.solver = 'CG';
        p.nIter  = 50;
        p.tol    = 10^-10;
end

CalculateDependentParameters();

    function Initialize()
        % counting parameters
        p.Ne   = 1; 
        p.Nc   = 1;
        p.Nsl  = 1;
        p.Nd   = 1;
        p.d_range = [];
        p.Nm   = [];
        p.Nx_o = [];
        p.Ny_o = [];
        p.Nx_i = [];
        p.Ny_i = [];
        p.Ns   = [];
        
        % FO parameters
        p.FO_e.flag = 1;
        p.FO_d.flag = 1;
        p.FOV   = [];
        p.res_o = [];
        p.res_i = [];
        
        % CS parameters
        p.CS_e.flag = 0;
        p.CS_d.flag = 0;
        p.TE0 = 0;
        p.dTE = 0;
        p.TE  = [];
        p.met = [];
        p.fm  = 0;
        
        % B0 parameters
        p.B0_e.flag   = 0;
        p.B0_e.vector = 0;
        p.B0_e.val    = [];
        p.B0_e.file   = [];
        p.B0_d.flag   = 0;
        p.B0_d.val    = 0;
        p.B0_d.file   = [];
        
        % through-plane B0 parameters
        p.tp_e.flag   = 0;
        p.tp_e.vector = 0;
        p.tp_e.val    = [];
        
        % T2* parameters
        p.T2_e.flag   = 0;
        p.T2_e.vector = Inf;
        p.T2_e.val    = [];
        p.T2_e.same   = 1; % 0: different T2* per met, 1: same T2* per met
        p.T2_d.flag   = 0;
        p.T2_d.vector = [];
        p.T2_d.val    = Inf;
        p.T2_d.same   = 1;
        
        % files
        p.obj_file  = [];
        p.traj_file = [];
        p.B0_file   = [];
        p.save_file = ['./Results/',num2str(p_no),'.mat'];
        
        % noise
        p.std = 0;
        
        % regularization
        p.lambda.vector = 0;
        p.lambda.val = [];
        
        % solver
        p.solver = [];
        p.nIter  = [];
        p.tol    = [];
    end

    function CalculateDependentParameters()
        p.Nm   = length(p.fm);
        p.TE   = p.TE0+(0:p.Ne-1)*p.dTE;
        p.Ntp  = length(p.tp_e.vector);
        p.Nx_o = round(p.FOV(1)/p.res_o(1));
        p.Ny_o = round(p.FOV(2)/p.res_o(2));
        p.Nx_i = round(p.FOV(1)/p.res_i(1));
        p.Ny_i = round(p.FOV(2)/p.res_i(2));
        
        % if no object is loaded: use point object per metabolite
        if isempty(p.obj_file) 
            p.Nd = p.Nm;
        end
        
        % if no d_range is defined: calculate all dynamics
        if isempty(p.d_range)
            p.d_range = 1:p.Nd;
        end
        
        % encoding B0 map always has same amount of slices as object
        p.B0_e.Nsl = p.Nsl;

        % if same T2* is used for all met: add values to T2* vector
        if p.T2_e.same == 1
            tmp = p.T2_e.vector;
            p.T2_e.vector = zeros(length(tmp)*p.Nm,1);
            for i = 1:length(tmp)
                p.T2_e.vector(p.Nm*(i-1)+1:p.Nm*i) = repmat(tmp(i),[p.Nm,1]);
            end
        end
    end

end

