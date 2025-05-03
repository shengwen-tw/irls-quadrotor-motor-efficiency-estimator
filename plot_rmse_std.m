rmse_irls = [0.023358, 0.0237185, 0.02357875];
rmse_ekf  = [0.02247225, 0.027893, 0.02250075];
std_irls  = [0.0233525, 0.02008725, 0.0235735];
std_ekf   = [0.02246, 0.02005375, 0.022483];

case_labels = {'Degrade only','Abrupt fault only','Combined case'};
y_label     = 'Value';
save_name   = 'rmse_std_subplots';

figure();

% RMSE graph
nexttile;
vals_rmse = [rmse_irls(:), rmse_ekf(:)];
bar(vals_rmse, 1, 'grouped'); grid on; hold on;  % <<< bar width=1

set(gca,'XTick',1:3,'XTickLabel',case_labels,'FontSize',11);
ylabel(y_label,'FontSize',11);
title('RMSE','FontSize',12);
legend({'IRLS','EKF'},'Location','northwest', 'Orientation','horizontal');

% Mark values above the bars
ngroups = size(vals_rmse,1);
nbars   = size(vals_rmse,2);
groupwidth = 1.0; %
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j = 1:ngroups
        yv = vals_rmse(j,i);
        text(x(j), yv, sprintf('%.4f', yv), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end
end
ylim([0, max(vals_rmse(:))*1.25]);

% Standard deviation
nexttile;
vals_std = [std_irls(:), std_ekf(:)];
bar(vals_std, 1, 'grouped');
grid on;
hold on;

set(gca,'XTick',1:3,'XTickLabel',case_labels,'FontSize',11);
ylabel(y_label,'FontSize',11);
title('Standard Deviation','FontSize',12);

% Mark values above the bars
ngroups = size(vals_std,1);
nbars   = size(vals_std,2);
groupwidth = 1.0;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j = 1:ngroups
        yv = vals_std(j,i);
        text(x(j), yv, sprintf('%.4f', yv), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end
end
ylim([0, max(vals_std(:))*1.25]);

exportgraphics(gcf, [save_name '.png'], 'Resolution', 300);