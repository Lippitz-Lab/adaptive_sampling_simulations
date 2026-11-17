% sets variables needed for calculations

function [omega, mu, gamma] = make_molecule(p)

    arguments
        p (1,5) double {mustBePositive}
    end

    omega1 = p(1);
    omega2 = p(2);
    mu1 = p(3);
    mu2 = p(4);
    gamma0 = p(5);

    omega = NaN(3);
    omega(1,2) = omega1;
    omega(1,3) = omega2;
    omega(2,3) = omega2 - omega1;
    
    mu = NaN(3);
    mu(1,2) = mu1;
    mu(1,3) = mu2; 
    
    % set lifetime to 1 ns;
    tau = 1e6; % 1 ns in fs
    gamma = gamma0 * ones(3);
    gamma(1:3+1:end) = 1 / tau;  % diagonal
end