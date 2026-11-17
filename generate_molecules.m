% simulate molecular 1d and 2d spectra

function generate_molecules(Nmol)

    arguments
        Nmol (1,1) double {mustBeInteger, mustBePositive}   %Number of molecules generated
    end

    rng(3);

    
    % rotating wave frequencies
    omega21 = 2.3; % rad/fs, approx. 1.5 eV
    omega43 = 2.3; % rad/fs, approx. 1.5 eV
    t2 = 10;       % fs

    np = 20;       % number of points in each dimension
    dt = 20;       % time step in fs
    t_list = (0.5 + (0:np-1)) .* dt;
    [t3, t1] = meshgrid(t_list, t_list);

    % omega1, omega2, mu1, mu2, gamma
    p0 = [2.35, 2.25, 1.0, 1.3, 0.03];       % mean
    pnoise = [0.02, 0.02, 0.2, 0.2, 0.005];  % std. dev.

    gamma0 = p0(end);      % mean of dephasing
    omega1 = 2:0.001:3;    % linear spectrum freq. axis

    % Preallocate output arrays
    spec_lin = zeros(Nmol, numel(omega1));
    trace_nr = complex(zeros(Nmol, np, np));
    trace_r = complex(zeros(Nmol, np, np));

    % iterate over molecules
    for mol = 1:Nmol

        % draw random molecular parameter
        pmol = p0 + pnoise .* randn(size(p0));      %randomize transition frequencies etc

        % generate required variables
        [omega, mu, gamma] = make_molecule(pmol);       %3x3 Matricies with energy levels,transition dipole moments, linewidths

        % calculate 
        [one_trace_nr, one_trace_r, one_spec_1d] = core_formulas( ...
            mu, omega, omega21, omega43, t1, t2, t3, gamma, omega1);    %Calculate rephasing, non rephasing (time domain traces) and linear spec from the model molecule

        % store data
        spec_lin(mol, :) = -real(one_spec_1d(:)).';

        % scale lowest coeff. of rephasing trace to abs. value of one
        scale = abs(one_trace_r(1, 1));
        trace_nr(mol, :, :) = one_trace_nr ./ scale;
        trace_r(mol, :, :) = one_trace_r ./ scale;
    end


    % save to file
    omega_lin = omega1;

    save("molecules.mat", ...
        "trace_r", "trace_nr", "t1", "t3", "t_list", "Nmol", ...
        "omega21", "omega43", "dt", "p0", "pnoise", ...
        "gamma0", "omega_lin", "spec_lin");

end