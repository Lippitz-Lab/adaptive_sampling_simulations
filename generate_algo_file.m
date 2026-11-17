% simulate the adaptive algorithm
% Here we use random numbers as shot noise of individual mesurements
% this takes some time! About 0.4 sec per molecule and noise level !

function generate_algo_file

    rng(2);

    disp("generate_algo_file");

    load("molecules.mat", "t1", "t3", "gamma0", "Nmol", "trace_r");
    load("daemon_and_square.mat", "noiselevels", "allchis", "Nmax");

    nNoise = numel(noiselevels);

    % Preallocate output arrays
    algo_chi = NaN(Nmol, nNoise, Nmax);
    algo_m = NaN(Nmol, nNoise, Nmax);

    for mol = 1: Nmol
    
        if mod(mol, 10) == 0
            disp(mol);
        end
        signal_truth = squeeze(trace_r(mol, :, :));

        for id = 1:nNoise
            [chilist, mlist, param_trace] = adaptive_algo_scan( ...
                signal_truth, noiselevels(id), Nmax, gamma0, ...
                t1 + t3, 0);

            algo_chi(mol, id, :) = chilist;
            algo_m(mol, id, :) = mlist;
            algo_amp(mol, id, :) = param_trace(:,1);
            algo_gamma(mol, id, :) = param_trace(:,2);
            algo_ym(mol, id, :) = param_trace(:,3);
            algo_yd(mol, id, :) = param_trace(:,4);
        end
    end


    algo_std = squeeze(std(algo_chi, 1, "omitnan"));
    [~, algo_std] = merge_noiselevels(algo_std, noiselevels);

    algo_chi = squeeze(mean(algo_chi, 1, "omitnan"));
    [~, algo_chi] = merge_noiselevels(algo_chi, noiselevels);

    algo_amp = squeeze(mean(algo_amp, 1, "omitnan"));
    [~, algo_amp] = merge_noiselevels(algo_amp, noiselevels);

    algo_ym = squeeze(mean(algo_ym, 1, "omitnan"));
    [~, algo_ym] = merge_noiselevels(algo_ym, noiselevels);
    
    algo_yd = squeeze(mean(algo_yd, 1, "omitnan"));
    [~, algo_yd] = merge_noiselevels(algo_yd, noiselevels);

    algo_gamma = squeeze(mean(algo_gamma, 1, "omitnan"));
    [~, algo_gamma] = merge_noiselevels(algo_gamma, noiselevels);

    algo_m = squeeze(mean(algo_m, 1, "omitnan"));
    [~, algo_m] = merge_noiselevels(algo_m, noiselevels);

    allchis.algo = algo_chi;
    allchis.algo_std = algo_std;
    allchis.algo_amp = algo_amp;
    allchis.algo_gamma = algo_gamma;
    allchis.algo_m = algo_m;

    allchis.algo_ym = algo_ym;
    allchis.algo_yd = algo_yd;

    
    save("algo_chi.mat", "allchis");

end