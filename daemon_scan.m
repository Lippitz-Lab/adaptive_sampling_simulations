% ideal Daemon

function [chilist, mlist, snapshot] = daemon_scan( ...
        signal_truth, noiselevel, Nmax, chi_snapshot)

    cn2 = abs(signal_truth(:)).^2;   % |c_n|^2
    signalNorm = sum(cn2, "all");
    noiseVariance = noiselevel^2;

    Nmeas = zeros(size(cn2));        % N_n
    snapshot = [];

    chilist = zeros(Nmax, 1);
    mlist = zeros(Nmax, 1);

    for N = 1:Nmax

        % reduction in chisq when measured at this index
        delta = cn2 ./ (cn2 .* Nmeas + noiseVariance) ...
            - cn2 ./ (cn2 .* (Nmeas + 1) + noiseVariance);

        % find maximum chis reduction and measure there
        [~, idnext] = max(delta);
        Nmeas(idnext) = Nmeas(idnext) + 1;

        %  optimal window function
        window = cn2 .* Nmeas ./ ...
            (cn2 .* Nmeas + noiseVariance);

        %  expectation value of residuals
        res = cn2;
        measured = Nmeas > 0;
        res(measured) = ...
            window(measured) .* noiseVariance ./ Nmeas(measured);

        % expectation value of MSE
        chisq = sum(res, "all") ./ signalNorm;

        % store progress 
        chilist(N) = chisq;
        mlist(N) = sum(measured, "all");  % number of unique coeffs.

        % store snapshot once, as soon as chi_snapshot is passed
        if isempty(snapshot) && chisq < chi_snapshot
            snapshot.Nmeas = reshape(Nmeas, size(signal_truth));
            snapshot.window = reshape(window, size(signal_truth));
            snapshot.residuals_sq = reshape( ...
                res ./ signalNorm, size(signal_truth));
            snapshot.chisq = chisq;
            snapshot.N = sum(Nmeas, "all");
        end
    end

    % store snapshot at the end if chi_snapshot == 0
    if  chi_snapshot == 0
        snapshot.Nmeas = reshape(Nmeas, size(signal_truth));
        snapshot.window = reshape(window, size(signal_truth));
        snapshot.residuals_sq = reshape( ...
            res ./ signalNorm, size(signal_truth));
        snapshot.chisq = chisq;
        snapshot.N = sum(Nmeas, "all");
    end

end