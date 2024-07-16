function highlightInputs(x, inputs, color, alpha, axis)

if size(color,1)==1
    color=[1,1,1; color];
end

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


%     if ~isempty(ons) && ~isempty(offs)
%         if max(ons) > max(offs)
%             offs = [offs, find(~isnan(inputs), 1, 'last')];
%         end
%         if min(offs) < min(ons)
%             ons = [find(~isnan(inputs), 1, 'first'), ons];
%         end
%     elseif isempty(ons) && ~isempty(offs)
%         ons = find(~isnan(inputs), 1, 'first');
%     elseif isempty(offs) && ~isempty(ons)
%         offs = find(~isnan(inputs), 1, 'last');
%     elseif (all(input_differences==0)&&inputs(1)>0)
%         ons=1;
%         offs=length(x);
%     else
%         return;
%     end
%
%     for i = 1:length(ons) - length(offs)
%         offs = [offs, find(~isnan(inputs), 1, 'last')];
%     end
%
%     for i = 1:length(ons)
%         x_start = x(ons(i));
%         x_end = x(offs(i));
%         a = area([x_start, x_end], [1,1]*axis.YLim(2));
%         a.FaceColor = color;
%         a.FaceAlpha = alpha*inputs(ons(i))/max(inputs);
%         a.EdgeAlpha = 0;
%         a.Tag = 'highlighted_area';
%     end
end

