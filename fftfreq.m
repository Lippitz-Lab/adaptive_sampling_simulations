function omega = fftfreq(dt, nsample)
% calculates omega axis from sampling time interval
% (and vice versa, just use dt as domega)

    n_half = fix(nsample/2); % round to integer
    n_remain = nsample - 2 * n_half;
    fs = abs(1. / dt); % sampling frequency
    freq = (-1 * n_half :1:n_half -1 + n_remain)./ nsample .* fs;
    omega = freq .* 2 .* pi;
    
    
end
    