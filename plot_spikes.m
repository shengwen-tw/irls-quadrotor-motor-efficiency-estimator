maxspike_irls_fault   = [0.139575, 0.123851, 0.142402, 0.130006];  % Abrupt fault only (IRLS)
maxspike_ekf_fault    = [0.286273, 0.263499, 0.292192, 0.274747];  % Abrupt fault only (EKF)

maxspike_irls_combo   = [0.120609, 0.102525, 0.129096, 0.150433];  % Degrade + Abrupt fault (IRLS)
maxspike_ekf_combo    = [0.176225, 0.155338, 0.259303, 0.220348];  % Degrade + Abrupt fault (EKF)

motor_labels = {'Motor 1','Motor 2','Motor 3','Motor 4'};
y_label      = '|Max Spike|';
save_name    = 'max_spike_per_motor';

figure();

% Abrupt fault only
ax1 = nexttile;
vals_fault = [maxspike_irls_fault(:) maxspike_ekf_fault(:)];   % 4x2
bar(vals_fault, 1, 'grouped');
grid(ax1,'on');
hold(ax1,'on');

set(ax1,'XTick',1:4,'XTickLabel',motor_labels,'FontSize',11);
ylabel(ax1, y_label, 'FontSize',11);
title(ax1, 'Abrupt fault only', 'FontSize',12);
lgd = legend(ax1, {'IRLS','EKF'}, 'Location','northwest', 'Orientation','horizontal');

ylim(ax1, [0, 0.43]);

% Mark values above the bars
ng = size(vals_fault,1); nb = size(vals_fault,2); gw = 1.0;
for i = 1:nb
    x = (1:ng) - gw/2 + (2*i-1)*gw/(2*nb);
    for j = 1:ng
        text(ax1, x(j), vals_fault(j,i), sprintf('%.4f', vals_fault(j,i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end
end

% Degrade + Abrupt fault
ax2 = nexttile;
vals_combo = [maxspike_irls_combo(:) maxspike_ekf_combo(:)];   % 4x2
bar(vals_combo, 1, 'grouped'); grid(ax2,'on'); hold(ax2,'on');

set(ax2,'XTick',1:4,'XTickLabel',motor_labels,'FontSize',11);
ylabel(ax2, y_label, 'FontSize',11);
title(ax2, 'Degrade and abrupt fault', 'FontSize',12);

% Mark values above the bars
ng = size(vals_combo,1); nb = size(vals_combo,2); gw = 1.0;
for i = 1:nb
    x = (1:ng) - gw/2 + (2*i-1)*gw/(2*nb);
    for j = 1:ng
        text(ax2, x(j), vals_combo(j,i), sprintf('%.4f', vals_combo(j,i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end
end

ymax_auto = max(vals_combo(:)) * 1.15;
ylim(ax2, [0, ymax_auto]);
exportgraphics(gcf, [save_name '.png'], 'Resolution', 300);