function params = nj_default_parameters(~)
%nj_default_parameters. Loads default parameters for Nora's Neuropixels data
%
% params = nj_default_parameters()
%
% 2024-2025, Alexander Heimel

%% General

params.projectfolder = 'C:\Users\heimel.HERSENINSTITUUT\OneDrive\Projects\Jamann_et_al';

%% Analysis

params.use_fixed_random_seed = true;

params.use_only_good_clusters = true;

params.stringent_post_area = true; % if true, only selects POm, otherwise all thalamic regions

params.stringent_cell_quality = true; % if true, only select bombcell GOOD quality clusters

params.dt = 0.002; % s, bin width

params.max_isi = 0.010; % s, max. time between spikes to be called a burst
   % was 0.020 but changed because 100 Hz was used in paper

params.bursting_threshold = 10; 
% if < 1, then min. burst fraction to be called bursting
% if >=1, then min. number of bursts to be called bursting

params.fast_burst_max_isi = 0.004; % >250 Hz  

params.compute_mean_burst_isi_as_inverse_freq = false; % 

params.analyze_pairs = true; % analyze possible connected pairs

params.max_delay = 0.020; % s, max. allowed difference between pre and post

params.n_jitters = 50; % number of jitters of spike cross-correlogram

params.min_corr_peak_in_std = 4; % min. height of spike cross correlogram peak rel. to jittered

params.isi_hist_edges = [0 logspace(-3,3,19)];

params.only_pair_on_first_spikes = false; % only use first spikes in burst of pre-synaptic neuron for pair determination

params.compute_individual_burst_success = true; % slow additional analysis

params.remove_spike_doublets = true; % removes spikes with identical spiketimes in each cluster

%% Results

params.show_pairs = true;

%% Figures
params.clr_ctrl = [0 0 0];
params.clr_cpz = [0 0 0 ];
params.clr_l5 = [0.94 0.35 0.16];
params.clr_pom = [0.14 0.45 0.70];

%% Load processparams_local. Keep at the end
if exist('processparams_local.m','file')
    params = processparams_local( params );
end
