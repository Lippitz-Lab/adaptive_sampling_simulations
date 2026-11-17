% generates chisq(N) curves for square patter and Daemon
%
% We take shotnoise into account by the expectation value of the noise
% distributions. We do not use random numbers here.
% Only the molecules are drawn from their distribution


function generate_square_and_daemon_files(noiselevels)

    arguments
        noiselevels double {mustBeFinite, mustBeNonnegative, mustBeVector}  %
    end

    noiselevels = noiselevels(:);

    disp("generate_square_and_daemon_files")

    % Load  the required variables
    load("molecules.mat", "trace_r", "Nmol");

    Nmax = 1e3;
    nNoise = numel(noiselevels);

    % Preallocate output arrays
    daemon_chi = NaN(Nmol, nNoise, Nmax);
    daemon_m = NaN(Nmol, nNoise, Nmax);
    square3 = NaN(Nmol, nNoise, Nmax);
    square5 = NaN(Nmol, nNoise, Nmax);
    square7 = NaN(Nmol, nNoise, Nmax);

    for mol = 1:Nmol
        signal_truth = abs(squeeze(trace_r(mol, :, :)));

        for id = 1:nNoise

            [chilist, mlist] = daemon_scan( ...
                signal_truth, noiselevels(id), Nmax, eps);

            square3list = square_scan( ...
                3, signal_truth, noiselevels(id), Nmax);
            square5list = square_scan( ...
                5, signal_truth, noiselevels(id), Nmax);
            square7list = square_scan( ...
                7, signal_truth, noiselevels(id), Nmax);

            daemon_chi(mol, id, :) = reshape(chilist, 1, 1, Nmax);
            daemon_m(mol, id, :) = reshape(mlist, 1, 1, Nmax);
            square3(mol, id, :) = reshape(square3list, 1, 1, Nmax);
            square5(mol, id, :) = reshape(square5list, 1, 1, Nmax);
            square7(mol, id, :) = reshape(square7list, 1, 1, Nmax);
        end
    end

    % Combine simulations of different noise levels and average ensemble

    square3_std = squeeze(std(square3, 1, "omitnan"));
    [~, square3_std] = merge_noiselevels(square3_std, noiselevels);

    square5_std = squeeze(std(square5, 1, "omitnan"));
    [~, square5_std] = merge_noiselevels(square5_std, noiselevels);

    square7_std = squeeze(std(square7, 1, "omitnan"));
    [~, square7_std] = merge_noiselevels(square7_std, noiselevels);

    square3 = squeeze(mean(square3, 1, "omitnan"));
    square5 = squeeze(mean(square5, 1, "omitnan"));
    square7 = squeeze(mean(square7, 1, "omitnan"));

    [~, square3] = merge_noiselevels(square3, noiselevels);
    [~, square5] = merge_noiselevels(square5, noiselevels);
    [~, square7] = merge_noiselevels(square7, noiselevels);

    allchis.square3 = square3;
    allchis.square3_std = square3_std;

    allchis.square5 = square5;
    allchis.square5_std = square5_std;
    
    allchis.square7 = square7;
    allchis.square7_std = square7_std;

    daemon_std = squeeze(std(daemon_chi, 1, "omitnan"));
    [~, daemon_std] = merge_noiselevels(daemon_std, noiselevels);

    daemon_chi = squeeze(mean(daemon_chi, 1, "omitnan"));
    [x, daemon_chi] = merge_noiselevels(daemon_chi, noiselevels);

    allchis.daemon = daemon_chi;
    allchis.daemon_std = daemon_std;
    allchis.x = x;

    daemon_m = squeeze(mean(daemon_m, 1, "omitnan"));
    [~, daemon_m] = merge_noiselevels(daemon_m, noiselevels);
    allchis.daemon_m = daemon_m;

    save("daemon_and_square.mat", "noiselevels", "allchis", "Nmax");

end