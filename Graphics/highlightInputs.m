function highlightInputs(x, inputs, color, alpha, axis)


    hold on
    inputs = fillmissing(inputs, 'constant', NaN);
    if size(inputs,1) > size(inputs,2); inputs=inputs'; end
    input_differences = [0,diff(inputs)];
    input_differences(isnan(input_differences)) = 0;

    if isempty(x)
        x = linspace(0, length(inputs)-1, length(inputs));
    end

    if nargin < 5 || isempty(axis)
        axis = gca;
    end

    on_value  = max(quantile(input_differences,0.75),0.001);
    off_value = min(quantile(input_differences,0.25),-0.001);

    ons = find(input_differences >= on_value);
    offs = find(input_differences <= off_value);

    if ~isempty(ons) && ~isempty(offs)
        if max(ons) > max(offs)
            offs = [offs, find(~isnan(inputs), 1, 'last')];
        end
        if min(offs) < min(ons)
            ons = [find(~isnan(inputs), 1, 'first'), ons];
        end
    elseif isempty(ons) && ~isempty(offs)
        ons = find(~isnan(inputs), 1, 'first');
    elseif isempty(offs) && ~isempty(ons)
        offs = find(~isnan(inputs), 1, 'last');
    else
        return;
    end

    for i = 1:length(ons) - length(offs)
        offs = [offs, find(~isnan(inputs), 1, 'last')];
    end

    ons = sort(ons);
    offs = sort(offs);

    for i = 1:length(ons)
        x_start = x(ons(i));
        x_end = x(offs(i));
        a = area([x_start, x_end], [1,1]*axis.YLim(2));
        a.FaceColor = color;
        a.FaceAlpha = alpha*inputs(ons(i))/max(inputs);
        a.EdgeAlpha = 0;
        a.Tag = 'highlighted_area';
    end
end

