

% posterior sigmas plot

% load all 26 outputs
load output_mat_files/all_combos/est01.mat;
load output_mat_files/all_combos/est02.mat;
load output_mat_files/all_combos/est03.mat;
load output_mat_files/all_combos/est04.mat;
load output_mat_files/all_combos/est05.mat;
load output_mat_files/all_combos/est06.mat;
load output_mat_files/all_combos/est07.mat;
load output_mat_files/all_combos/est08.mat;
load output_mat_files/all_combos/est09.mat;
load output_mat_files/all_combos/est10.mat;

load output_mat_files/all_combos/est11.mat;
load output_mat_files/all_combos/est12.mat;
load output_mat_files/all_combos/est13.mat;
load output_mat_files/all_combos/est14.mat;
load output_mat_files/all_combos/est15.mat;
load output_mat_files/all_combos/est16.mat;
load output_mat_files/all_combos/est17.mat;
load output_mat_files/all_combos/est18.mat;
load output_mat_files/all_combos/est19.mat;
load output_mat_files/all_combos/est20.mat;

load output_mat_files/all_combos/est21.mat;
load output_mat_files/all_combos/est22.mat;
load output_mat_files/all_combos/est23.mat;
load output_mat_files/all_combos/est24.mat;
load output_mat_files/all_combos/est25.mat;
load output_mat_files/all_combos/est26.mat;

nD = length(est26);

for i = 1:nD
    input(1).est(i) = est01(i);
    input(2).est(i) = est02(i);
    input(3).est(i) = est03(i);
    input(4).est(i) = est04(i);
    input(5).est(i) = est05(i);
    input(6).est(i) = est06(i);
    input(7).est(i) = est07(i);
    input(8).est(i) = est08(i);
    input(9).est(i) = est09(i);
    input(10).est(i) = est10(i);
    input(11).est(i) = est11(i);
    input(12).est(i) = est12(i);
    input(13).est(i) = est13(i);
    input(14).est(i) = est14(i);
    input(15).est(i) = est15(i);
    input(16).est(i) = est16(i);
    input(17).est(i) = est17(i);
    input(18).est(i) = est18(i);
    input(19).est(i) = est19(i);
    input(20).est(i) = est20(i);
    input(21).est(i) = est21(i);
    input(22).est(i) = est22(i);
    input(23).est(i) = est23(i);
    input(24).est(i) = est24(i);
    input(25).est(i) = est25(i);
    input(26).est(i) = est26(i);
end

for k = 1:26
    est = input(k).est;
    for i = 1:nD
        inpt(1).sig(k,i) = est(i).eTC;
        inpt(2).sig(k,i) = est(i).eTA;
        inpt(3).sig(k,i) = est(i).tp(1).eph;
        inpt(4).sig(k,i) = est(i).tp(2).epco2;
        inpt(5).sig(k,i) = est(i).tp(3).eco3;
    end
end
%%
% reorder to same order as f_hat plot
for i = 1:5
    sig = inpt(i).sig;
    sig_ord(1,:) = sig(9,:);
    sig_ord(2,:) = sig(8,:);
    sig_ord(3,:) = sig(20,:);
    sig_ord(4,:) = sig(10,:);
    sig_ord(5,:) = sig(7,:);
    sig_ord(6,:) = sig(16,:);
    sig_ord(7,:) = sig(3,:);
    sig_ord(8,:) = sig(1,:);
    sig_ord(9,:) = sig(19,:);
    sig_ord(10,:) = sig(13,:);
    sig_ord(11,:) = sig(6,:);
    sig_ord(12,:) = sig(12,:);
    sig_ord(13,:) = sig(23,:);
    sig_ord(14,:) = sig(15,:);
    sig_ord(15,:) = sig(4,:);
    sig_ord(16,:) = sig(2,:);
    sig_ord(17,:) = sig(18,:);
    sig_ord(18,:) = sig(5,:);
    sig_ord(19,:) = sig(24,:);
    sig_ord(20,:) = sig(11,:);
    sig_ord(21,:) = sig(14,:);
    sig_ord(22,:) = sig(22,:);
    sig_ord(23,:) = sig(25,:);
    sig_ord(24,:) = sig(17,:);
    sig_ord(25,:) = sig(21,:);
    sig_ord(26,:) = sig(26,:);
    outpt(i).sig_ord = sig_ord;
end


% next, plot the posterior sigmas

lbl = {'pH, $p$CO$_{2}$',...            % 9 
    'pH, CO$_{3}$',...                  % 8 
    'pH, $p$CO$_{2}$, CO$_{3}$',...     % 20
    '$p$CO$_{2}$, CO$_{3}$',...         % 10
    'A$_T$, CO$_{3}$',...               % 7
    'C$_T$, $p$CO$_{2}$, CO$_{3}$',...  % 16
    'C$_T$, $p$CO$_{2}$',...            % 3 
    'C$_T$, A$_T$',...                  % 1 
    'A$_T$, $p$CO$_{2}$, CO$_{3}$',...  % 19 
    'C$_T$, A$_T$, CO$_{3}$',...        % 13
    'A$_T$, $p$CO$_{2}$',...            % 6 
    'C$_T$, A$_T$, $p$CO$_{2}$',...     % 12
    'C$_T$, A$_T$, $p$CO$_{2}$, CO$_{3}$',...   % 23
    'C$_T$, pH, CO$_{3}$',...                   % 15
    'C$_T$, CO$_{3}$', ...                      % 4 
    'C$_T$, pH',...                             % 2
    'A$_T$, pH, CO$_{3}$',...                   % 18
    'A$_T$, pH',...                             % 5
    'C$_T$, pH, $p$CO$_{2}$, CO$_{3}$',...      % 24
    'C$_T$, A$_T$, pH',...                      % 11 
    'C$_T$, pH, $p$CO$_{2}$',...                % 14
    'C$_T$, A$_T$, pH, CO$_{3}$',...            % 22
    'A$_T$, pH, $p$CO$_{2}$, CO$_{3}$',...      % 25
    'A$_T$, pH, $p$CO$_{2}$',...                % 17   
    'C$_T$, A$_T$, pH, $p$CO$_{2}$',...         % 21
    'C$_T$, A$_T$, pH, $p$CO$_{2}$, CO$_{3}$'}; % 26

% AXES PROPERTIES
set(groot,'defaultAxesFontName','Athelas',...
    'defaultAxesFontSize',8,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','off',...
    'defaultAxesYMinorTick','off');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Athelas',...
    'defaultTextInterpreter','latex');

ddgreen = [0, 102/255, 0];
clr = ddgreen;

%% 

t = tiledlayout(5,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
b = boxchart((outpt(1).sig_ord)');
b.JitterOutliers = 'off';
b.MarkerStyle = '.';
b.MarkerSize = 4;
b.BoxWidth = 0.6;
b.LineWidth = 1.25;
b.WhiskerLineColor = clr;
b.BoxFaceColor = clr;
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
ylabel('$\mathbf{\sigma_{post}}$ for C$_T$') % normal
ax = gca;
ax.FontSize = 7;
grid on;

nexttile
b = boxchart((outpt(2).sig_ord)');
b.JitterOutliers = 'off';
b.MarkerStyle = '.';
b.MarkerSize = 4;
b.BoxWidth = 0.6;
b.LineWidth = 1.25;
b.WhiskerLineColor = clr;
b.BoxFaceColor = clr;
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
ylabel('$\mathbf{\sigma_{post}}$ for A$_T$') % normal
ax = gca;
ax.FontSize = 7;
grid on;

nexttile
b = boxchart((outpt(3).sig_ord)');
b.JitterOutliers = 'off';
b.MarkerStyle = '.';
b.MarkerSize = 4;
b.BoxWidth = 0.6;
b.LineWidth = 1.25;
b.WhiskerLineColor = clr;
b.BoxFaceColor = clr;
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
ylabel('$\mathbf{\sigma_{post}}$ for $p$H')
ax = gca;
ax.FontSize = 6;
grid on;

nexttile
b = boxchart((outpt(4).sig_ord)');
b.JitterOutliers = 'off';
b.MarkerStyle = '.';
b.MarkerSize = 4;
b.BoxWidth = 0.6;
b.LineWidth = 1.25;
b.WhiskerLineColor = clr;
b.BoxFaceColor = clr;
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
ylabel('$\mathbf{\sigma_{post}}$ for $p$CO$_2$') % normal
ax = gca;
ax.FontSize = 7;
grid on;

nexttile
b = boxchart((outpt(5).sig_ord)');
b.JitterOutliers = 'off';
b.MarkerStyle = '.';
b.MarkerSize = 4;
b.BoxWidth = 0.6;
b.LineWidth = 1.25;
b.WhiskerLineColor = clr;
b.BoxFaceColor = clr;
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
ylabel('$\mathbf{\sigma_{post}}$ for CO$_3^{2-}$','FontSize',7) % normal
ax = gca;
ax.FontSize = 6.5;
grid on;
xticklabels(lbl);



h = gcf; 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);

print(h,'posterior_sigmas_plot.pdf','-dpdf');






