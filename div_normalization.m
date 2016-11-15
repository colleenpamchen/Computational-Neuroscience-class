function div_normalization()

    % Script for simulating divisive normalization in visual cortex
    % corresponding to orientation
    %
    % There are seven neurons being considered corresponding to the
    % orientation which is most likely to make them fire the most given a
    % stimulus in the direction as [0, 22.5, 45, 90, 112.5, 135, 180] in
    % degrees.
    % These neurons receive inputs from seven corresponding input neurons
    % as well as one inhibitory neuron
    % The inhibitory neuron in turn is connected to all the input neurons
    % The presence of inhibitory neuron gives the effect of divisive
    % normalization as given by the equation 8 in Carandini & Heeger
    % 
    % I use Gaussian tuning curves with maximum response as 52.14 Hz and
    % standard deviation to be 14.73 degrees
    %
    % There are two stimuli: grating with 0 degree orientation and grating
    % with 90 degree orientation with varying percentages of contrast
    % from 0% to 100%
    %
    % I plot a 5x5 matrix where each row corresponds to a contrast
    % percentage for grating with orientation 0 and each column with a
    % contrast percentage for grating with orientation 90
    %
    % The plots show that when only one of the two stimuli are present, the
    % neuron sensitive to the particular orientation has a higher mean
    % firing rate. When both are present, they serve as noise for each
    % other and bring down the firing rate for the neuron sensitive to the
    % orientation with higher contrast. This symbolizes the effect of
    % normalization, where an oriented stimulus that is not supposed to
    % contribute to the firing rate for a particular neuron, results in 
    % reducing its firing rate if present as in sufficiently high contrast
    %
    % I plot another 5x5 matrix with similar plots in each cell, the
    % difference being that I take out the normalization effect of the
    % inhibitory neuron from the network and as expected, it can be seen
    % that the presence of an orthogonal stimulus does not affect the mean
    % firing rate of any neuron.

    r_max_gauss = 52.14;
    sigma_gauss = 14.73;
    function response = gaussian_tuning(s, in)
        response = r_max_gauss * exp(-0.5 * (in - s)^2 / sigma_gauss^2);
    end

    input = 90;
    cross_input = 0;
    
    neurons = [0, 22.5, 45, 90, 112.5, 135, 180];
    
    gaussian_responses = zeros(1, length(neurons));
    gaussian_cross_responses = zeros(1, length(neurons));
    
    for i = 1:length(neurons)
        s = neurons(i);
        gaussian_responses(i) = gaussian_tuning(s, input);
        gaussian_cross_responses(i) = gaussian_tuning(s, cross_input);
    end

    contrasts = [0, 25, 50, 75, 100];
    cross_contrasts = contrasts;
    
    gaussian = cell(length(cross_contrasts), length(contrasts));
    gaussian_normalized = cell(length(cross_contrasts), length(contrasts));
    
    sigma = 10;
    n = 2;
    
    for i = 1:length(contrasts)
        for j = 1:length(cross_contrasts)
            c_t = contrasts(i);
            c_m = cross_contrasts(j);
            norm = sigma^n + c_t^n + c_m^n;
            
            gaussian{j, i} = gaussian_responses * c_t^n + gaussian_cross_responses * c_m^n;
            gaussian_normalized{j, i} = gaussian{j, i} ./ norm;
            gaussian{j, i} = gaussian{j, i} ./ sigma^n;
        end
    end
    
    function plot_activity(x, y_cell, contrasts, cross_contrasts, name, y_limits)
        n_rows = length(cross_contrasts);
        n_cols = length(contrasts);
        
        figure();
        xq = 0:10:180;
        for r = 1:n_rows
            for c = 1:n_cols
                subplot(n_rows, n_cols, (r-1)*n_cols + c);
                y = y_cell{r, c};
                yq = interp1(x, y, xq, 'spline');
                plot(x, y, '.r', xq, yq, '-k');
                if exist('y_limits', 'var')
                    ylim(y_limits);
                end
                set(gca, 'XTick', [0, 90]);
                set(gca, 'XTickLabel', [0, 90]);
                if r == 1
                    title(strcat(int2str(contrasts(c)), '%'));
                end
                if c == 1
                    ylabel(strcat(int2str(cross_contrasts(r)), '%'), 'FontWeight','bold');
                end
            end
        end
        set(gcf,'name',name,'numbertitle','off')
    end

    plot_activity(neurons, gaussian, contrasts, cross_contrasts, 'Gaussian Tuning without Normalization', [-300, 5500]);
    plot_activity(neurons, gaussian_normalized, contrasts, cross_contrasts, 'Gaussian Tuning Normalized', [-10, 55]);
end

