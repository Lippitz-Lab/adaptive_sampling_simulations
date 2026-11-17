function [data_out, omega1, omega2] = my2Dfft_2axes_Trapezoidal(data_in, time1, time2, sf)
% 2D FFT with zero padding for t < 0
%        and omega axis

% we know that there is no signal (in this Feynman diagram) 
% for a different time order = negative times
% here we include this information by appending zeros

% time starts at zero for the first pixel and then increases
% t<0 pixel are thus in the second half of the array
%
% time1 is along first dimension of array = y-direction of imagesc =
% columns
    
    % trapezoidal rule (see Hamm Book)
    % each pixel contributes only with its fraction in quadrant 1
  
   [fac2xy , fac1xy ] = meshgrid(trap(time2) , trap(time1));
   data_in = data_in .* fac1xy .* fac2xy;

    % zero padding
    n1 = length(time1);
    n2 = length(time2);
    padded = zeros(sf*n1 + 1, sf*n2 + 1);
    padded(1:n1, 1:n2) = data_in ; %.* window ;
    
    data_out =fftshift(fft2(padded));
    
    omega1 = fftfreq( time1(2) - time1(1), sf*n1 + 1);
    omega2 = fftfreq( time2(2) - time2(1), sf*n2 + 1);

        
    % up to here, we assumed that the bottom-left element is at 
    % delay (0,0). Now we shift the spectrum to correct this
    
    [omega_2xy , omega_1xy ] = meshgrid(omega2 , omega1);
    d1 = time1(1);
    d2 = time2(1);
    
    data_out = data_out .* exp(-1i .* (d1 .*  omega_1xy + d2 .* omega_2xy));

end


function w = trap(time)
% make trapezoidal weights

   w = 0.5 +  time ./ (time(2) - time(1));
    
   [fm2, id2] = min(w);

    w(w > 1) = 1;
    w(w < 0) = 0;
    
    % first pixel too late
    if (fm2 > 1)
        w(id2) = fm2;
        disp('error: time axis starts too late, trying to correct');
    end
    
end