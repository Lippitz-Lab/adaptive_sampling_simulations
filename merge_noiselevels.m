% merges the different noiselevels together (to reduce simulation time)

function [x, y] = merge_noiselevels(data, noiselevels)

    data = squeeze(data);
    s = size(data);

    xcol = [];
    ycol = [];

    for id = 1:s(1)

        preclip = 101;
        
        x = (preclip:s(2));    
        
        xcol = [xcol, x ./ noiselevels(id)^2];
        ycol = [ycol, data(id,preclip:end)];
    end

    x = 10.^(1:0.01:8);
    y = interp1(xcol, ycol, x);

end
