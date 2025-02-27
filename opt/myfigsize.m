%% Plot commands
% figure;tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); nexttile;
% myfigstyle(gcf, 6, 6, 7, 7);
% saveas(gcf,'myfigure.fig');saveas(gcf,'myfigure.pdf');

%% Instruction
% > Measure real text width in cm
% > call myfigstyle with desired width, e.g. (12cm, 6cm)
% > for 2 plots besides: (6,6) and (6,6) with same font settings!

function myfigsize(fig, w, h, md, sm, lineWidth)
    
    set(fig,'Units','centimeters');
    xticks('manual');
    yticks('manual');
    pos = [0,0]; %get(fig,'Position');
    s = [0,0,0,0]; % additional padding

    set(fig,...
      'Position',[pos(1:2) w h],...
      'PaperUnits','centimeters',...
      'PaperPosition',[s(1)*w s(4)*h w h],...
      'PaperSize',[w*(1+s(1)+s(2)) h*(1+s(3)+s(4))]);
  
    if ~isnan(md) && ~isnan(sm)
        % Axes
        Axes = findall(fig, 'type', 'axes');
        numPlots = length(Axes);
        % if numPlots==1
        %     lineWidth = 1.5;
        %     markerSize = 6;
        %     markerSizeError = 12;      
        % else
        %     lineWidth = 2;
        %     markerSize = 4;
        %     markerSizeError = 15; % for dot '.', otherwise too large          
        % end

        % set(Axes, 'LineWidth', lineWidth/2);
        set(Axes, 'TickLabelInterpreter', 'latex');
        set(Axes, 'FontSize', md)

        Line = findall(fig, 'type', 'line');
        % for i=1:length(Line)
        %     if Line(i).LineWidth~=0
        %             Line(i).LineWidth=lineWidth;
        %     end
        %     if contains(Line(i).Marker,'.')
        %         if numPlots==1
        %             Line(i).MarkerSize = 2*markerSize; % double size for dot
        %         else
        %             Line(i).MarkerSize = 3*markerSize;
        %         end
        %     else
        %         Line(i).MarkerSize = markerSize;
        %     end
        % 
        % end
        if ~isnan(lineWidth)
            set(Line, 'LineWidth', lineWidth);
        end
        %set(Line, 'MarkerSize', markerSize);

        % Errorbar = findall(fig, 'type', 'errorbar');
        % set(Errorbar, 'LineWidth', lineWidth);
        % set(Errorbar, 'Marker', '.');
        % set(Errorbar, 'MarkerSize', markerSizeError);

        Intpr = findall(fig, '-property', 'Interpreter');
        set(Intpr, 'Interpreter', 'latex');

        Fonts = findall(fig, '-property', 'FontSize');
        set(Fonts, 'FontSize', md);

        % Legend Linewidth overwritte from findall('type', 'line')
        Legend = findall(fig, 'type', 'legend');
        set(Legend, 'FontSize', sm);
        % set(Legend, 'LineWidth', lineWidth/2);
        
        %% Colorbar
        cbar = findall(fig, 'type', 'colorbar');
        for c=cbar
            c.Label.Interpreter = 'Latex';
            c.TickLabelInterpreter = 'Latex';
            c.FontSize = md;
        end
        
    end

end
