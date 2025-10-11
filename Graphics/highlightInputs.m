function highlightInputs(x, inputs, color, alpha, axis)

if ischar(color)
    color = validatecolor(color);
end

if size(color,1)==1
    color=[1,1,1; color];
end

if ~strcmp(inputs,'None')
    hold on
    inputs = fillmissing(inputs, 'constant', NaN);
    
    if size(inputs,1) > size(inputs,2); inputs=inputs'; end
    
    inputs=round(inputs,2);
    
    input_differences = [0,diff(inputs)];
    input_differences(isnan(input_differences)) = 0;
    
    if isempty(x)
        x = linspace(0, length(inputs)-1, length(inputs));
    end
    
    if nargin < 5 || isempty(axis)
        axis = gca;
    end
    
    in_max=max(inputs);
    
    switches = find(abs(input_differences) > (in_max/100));
    switches = [1,switches,length(x)];
    
    if in_max>0
        for i = 1:length(switches)-1
            a = area([x(switches(i)), x(switches(i+1))], [1,1]*axis.YLim(2));
            norm_input = inputs(switches(i))/in_max;
            a.FaceColor = norm_input * color(end,:) + (1-norm_input) * color(1,:);
            a.FaceAlpha = alpha;
            a.EdgeAlpha = 0;
            a.Tag = 'highlighted_area';
        end
    end
end
end

