close all
clear all

% use export_fig from https://github.com/altmany/export_fig
% to export heatmaps as bitmaps with vector graphics labeling
% or else use the included exportgraphics
% to export heatmap and labels both as vector graphics
% (this gives problems in some pdf viewers)

use_export_fig = true;
if use_export_fig
    addpath("../export_fig/")
end

% intermediate data is stored as mat files. 
% comment out generate_XYZ calls to use existing data and save time.

Nmol = 1000;
generate_molecules(Nmol)    %Generates Nmol molecules and saves time domain trace of rephasing and non rephasing for each molecule, linear spectrum and other parameters like rotating frame frequencies, time steps, std for transition frequencies etc


% Generate the data for fig.3 (Daemon and square MSEs for the noiselevels;
% Merge the different noiselevels together (to reduce simulation time)
noiselevels =  sqrt(10.^(2:-1:-6)); 
generate_square_and_daemon_files(noiselevels) 
tic
generate_algo_file
Tmol = toc / (Nmol * length(noiselevels))


fig_chi_over_N % plots fig. 3 in the paper

% select demo molecule for fig. 2 and 4
idmol = 2;
% 
fig_mol_spectra(idmol, use_export_fig) % plots Fig. 2 in paper

generate_compare_snapshots_at_N(idmol)
fig_compare_both_at_N(idmol, use_export_fig) % plots Fig. 4 in paper


demo_fig1  % plots Fig. 1 in paper


