function create_task_screen(cue_pos, mov_dir, pre_dir, task_phase, varargin)
% create_task_screen Generate replica of task screen with extra info added
% Image generated is not to scale. Fixation/cue size and distance from
% center are not to scale.
% 
% INPUTS:
%   cue_pos:                        Angular position of cue
%   mov_dir:                        Angular position of movement
%   pre_dir:                        Angular position of predicted movement
%   task_phase:                     What task state are we currently in?
%   varargin
%       show_center_fixation:       Do we want to show the center fixation 
%                                   cue?
% OUTPUTS - None

%% Variable input parsing
p = inputParser;
p.addOptional('show_center_fixation', true);
p.parse(varargin{:});
inputs = p.Results;

%% Scaling so distance is correct on the plot.
% Since using `axis image`, the final plot is a square, but the main plot
% area ([-30 30 -30 30]) will be not a square.
cue_pos(1) = cue_pos(1) * 0.6;
mov_dir(1) = mov_dir(1) * 0.6;
pre_dir(1) = pre_dir(1) * 0.6;

%% Plot target cue
cla;
hold on;

if inputs.show_center_fixation
    plot(0, 0, 'rd', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 15);
end


if all(~isnan(cue_pos))
    plot(cue_pos(1), cue_pos(2), 'rd', 'MarkerFaceColor', 'r', ...
        'MarkerSize', 15);
end
%% Plot multiple things

if all(~isnan(mov_dir))
    quiver(0, 0, mov_dir(1), mov_dir(2), 'b', ...
            'MaxHeadSize', 5, ...
            'LineWidth', 5);
end

if all(~isnan(pre_dir))
   if all(~isnan(mov_dir)) && isequal(mov_dir, pre_dir)
       quiver(0, 0, pre_dir(1)*0.8, pre_dir(2)*0.8, 'w', ...
            'MaxHeadSize', 5, ...
            'LineWidth', 5);
   else
       quiver(0, 0, pre_dir(1), pre_dir(2), 'w', ...
            'MaxHeadSize', 5, ...
            'LineWidth', 5);
   end
end

hold off;
xlim([-30 30]);
ylim([-50 50]);
axis square;

% Add text
text(0, 40, task_phase, ...
    'FontSize', 18, ...
    'Color', 'w', ...
    'HorizontalAlignment', 'center');
text(0, -30, 'prediction', ...
    'FontSize', 15, ...
    'Color', 'w', ...
    'HorizontalAlignment', 'center');
text(0, -40, 'movement', ...
    'FontSize', 15, ...
    'Color', 'b', ...
    'HorizontalAlignment', 'center');

set(gca,'Color','k') % Background color
% Remove axes labels
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])


end
