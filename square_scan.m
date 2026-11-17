% simulate N x N square pattern

function chilist = square_scan(edgelength, signal_truth, noiselevel, Nmax)

    signal_truth = squeeze(signal_truth);

    % Generate square index pattern
    [x, y] = ndgrid(1:edgelength);
    p = [x(:), y(:)];

    p = repmat(p, ceil(Nmax / size(p, 1)), 1);
    id = sub2ind(size(signal_truth), p(:, 1), p(:, 2));

    % Initialize res_sq by |c_n|^2
    res_sq = abs(signal_truth).^2;
    signalNorm = sum(res_sq, "all");
    noiseVariance = noiselevel^2;

    chilist = zeros(Nmax, 1);
    Nmeas = zeros(size(signal_truth));

    for N = 1:Nmax
        Nmeas(id(N)) = Nmeas(id(N)) + 1;

        measured = Nmeas > 0;
        res_sq(measured) = noiseVariance ./ Nmeas(measured);

        chilist(N) = sum(res_sq, "all") ./ signalNorm;
    end

end