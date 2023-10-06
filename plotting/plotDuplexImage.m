function my_map = plotDuplexImage(X_img_mm, Z_img_mm, foreground, background,varargin)
% my_map = plotDuplexImage(X_img_mm, Z_img_mm, foreground, background,varargin)
%
% this function will plot the foreground as a mask on top of the
% background. It is a bit fragile so YMMV with the effectiveness of this
% function.
%
% Example usage:
% pixelsize = 0.1;
% X_img_mm = pixelsize/2 + (0:size(iDopP,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
% Z_img_mm = pixelsize/2 + (0:size(iDopP,1)-1)*pixelsize + UF.Depth(1);
% cmap = plotDuplexImage(X_img_mm, Z_img_mm, statistical_overlay, wallpaper);
%
%
% vargin:
% -----------------------------------------------------------
% X_img_mm: x axis indices in mm
% Z_img_mm: z axis indices in mm
% foreground: (zPix x xPix) array/image, note: NaNs will be 'transparent'
% background: (zPix x xPix) array/image
% 
% optional vargin:
% -----------------------------------------------------------
% colormap2use          (colormap array, default: hot)
% FontSize              (numeric, default: 18)
% flipMap               (bool, default: false) 
% nonlinear_fg          (bool, default: false), e.g. fg = nthroot(fg,4);
% nonlinear_bg          (bool, default: false), scaling for background
% scalebar              (numeric, default: 0 does not plot a scalebar, otherwise 
%                       numeric value for 1 mm, e.g. a value of 40 means 40px=1mm)
% showColorbar          (bool, default: false)
% AutoColorBarLimits    (bool/[min max], default: true). If not true, use [min
%                       max] to indicate colorbar limits.
% displayBest           (bool or double fraction, default: false. If you pass a
%                       double=1, it will print the bottom 10%. Otherwise if you
%                       pass a double fraction, e.g. 'displayBest',0.2, it will
%                       print that fraction, e.g. bottom 20% (usually p values)
% displayBot            same as displayBest
% displayTop            same as displayBest and displayBot except it will print the
%                       top quantile defined rather than the bottom.
% transparencyMap       (zPix x xPix) An image that is the same size as foreground 
%                       where 1 = opaque and 0 = transparent
% ColorBarTitle         (str) A title for the colorbar
% darkestShade          (int) Must be between [0 1). Adjusts darkest shade
%                       for the background image.
%
% optional return arguments:
% my_map: colormap

p = inputParser;
p.addParameter('colormap2use',hot());
p.addParameter('FontSize', 18, @isnumeric);
p.addParameter('flipMap',false);
p.addParameter('nonlinear_fg', false);
p.addParameter('nonlinear_bg', false);
p.addParameter('scalebar', 0, @isnumeric);
p.addParameter('showColorbar', false);
p.addParameter('AutoColorBarLimits',true);
p.addParameter('displayBest',false)
p.addParameter('displayBot',false)
p.addParameter('displayTop',false)
p.addParameter('transparencyMap', 1);
p.addParameter('ColorBarTitle', '');
p.addParameter('darkestShade', 0);
p.parse(varargin{:});
arg = p.Results;


% flip maps (if required)
if arg.flipMap
    foreground = flipud(foreground);
    background = flipud(background);
end

% normalize background -1:0
bgm = min(background(:));
bgM = max(background(:));
Ibg = (background - bgm)/(bgM-bgm)-1;

% normalize foreground 0:1
% If user provides colorbar limits, use those in place of automatic limits
if length(arg.AutoColorBarLimits)>1
    fgm = arg.AutoColorBarLimits(1);
    fgM = arg.AutoColorBarLimits(2);
    
else
    fgm = min(foreground(:));
    fgM = max(foreground(:));
end
Ifg = (foreground - fgm)/(fgM-fgm);


% If any value is greater than limits, set to max.
Ifg(Ifg>1) = 1;

% If any value is less than limits, set to min.
Ifg(Ifg<0) = 0;

% nonlinear scaling for better visualization
if arg.nonlinear_fg
    if arg.nonlinear_fg==1
        Ifg = nthroot(Ifg,5);
    else
        Ifg = nthroot(Ifg,arg.nonlinear_fg);
    end
end
if arg.nonlinear_bg
    if arg.nonlinear_bg==1
        Ibg = nthroot(Ibg+1,5);
    else
        Ibg = nthroot(Ibg+1,arg.nonlinear_bg);
    end
    Ibg = Ibg-1;
end

% displayBest, displayTop, displayBot
if arg.displayBot, arg.displayBest = arg.displayBot; end
% handle errors
if min([arg.displayBot, arg.displayBest, arg.displayTop]) < 0 ||...
        max([arg.displayBot, arg.displayBest, arg.displayTop]) > 1
    warning('display cutoff invalid'), return
end
% displaying both top & bottom if asked
if arg.displayBest && arg.displayTop 
    if arg.displayBest==1, arg.displayBest = 0.1; end
    if arg.displayTop==1, arg.displayTop = 0.1; end
    thresholdBest = quantile(Ifg(:),arg.displayBest);
    thresholdTop = quantile(Ifg(:),1-arg.displayTop);
    if isnan(thresholdBest), thresholdBest = fgm; end
    if isnan(thresholdTop),  thresholdTop  = fgM; end
    Ifg(Ifg>thresholdBest & Ifg<thresholdTop) = nan;    
elseif arg.displayBest 
    if arg.displayBest==1, arg.displayBest = 0.1; end    
    threshold = quantile(Ifg(:),arg.displayBest);    
    if ~isnan(threshold)
        Ifg(Ifg>threshold) = nan;         
    end    
elseif arg.displayTop
    if arg.displayTop==1, arg.displayTop = 0.1; end    
    threshold = quantile(Ifg(:),1-arg.displayTop);    
    if ~isnan(threshold)
        Ifg(Ifg<threshold) = nan;         
    end    
end
        

% combine the background & foreground
if all(size(Ifg)==size(Ibg))
    the_image = Ibg;
    the_image(~isnan(Ifg)) = Ifg(~isnan(Ifg));
    the_image(Ibg==min(Ibg(:))) = -1; % Sack a pixel to fix colormap
    %the_image(end, end) = -1; % Alternative option of pixel to sack.
    

    if any(~isnan(Ifg), 'all')
        % if any of the pixels in the foreground are not NaN, then the
        % colormap should go from -1:1.
        the_image(Ibg==max(Ibg(:))) = 1; % sac a pixel to fix colormap
        %the_image(1, 1) = 1; % Alternative option of a pixel to sack.
    else
        % If all pixels in the foreground are NaNs, then the colormap should go
        % from -1:0, not -1:1. 
        the_image(Ibg==max(Ibg(:))) = 0;
        %the_image(1, 1) = 0; % Alternative option of a pixel to sack.
    end
    
else
    error('the background & foreground sizes don''t match.')
end

% define colormap
base_map = arg.colormap2use;
my_map = ones(200,3);

if nnz(~isnan(Ifg)) ~= 0
    my_map(1:100,:) = repmat(linspace(arg.darkestShade,1,100)',[1 3]);         % b/w bg color
    my_map(101:200,:) = abs(min(imresize(base_map,[100,3]),1)); % fg colormap
else
    %Entire foreground is NaN
    my_map = repmat(linspace(arg.darkestShade,1,100)',[1 3]);
end


% plot result
% Handle transparency
if numel(arg.transparencyMap) == 1
    %If no transparency map, then don't play around with AlphaData to
    %preserve backwards compatibility
    imagesc(X_img_mm, Z_img_mm, the_image);
else
    % If a transparency map is provided, use it.
    % Overlay the transparent image on top of a solid image. This allows
    % the transparency to actually function correctly. Otherwise, the
    % background for high transparency areas is white.
    image_bg = Ibg;
    image_bg(Ibg==min(Ibg(:))) = -1; % sac a pixel to fix colormap
    if any(~isnan(Ifg), 'all')
        image_bg(Ibg==max(Ibg(:))) = 1;
    else
        image_bg(Ibg==max(Ibg(:))) = 0;
    end
    

    imagesc(X_img_mm, Z_img_mm, image_bg);
    % To lighten the background slightly
    %imagesc(X_img_mm, Z_img_mm, image_bg, 'AlphaData', 0.75 * ones(size(image_bg)));
    
    if any(~isnan(Ifg), 'all')
        hold on;
        image_fg = Ifg;

        imagesc(X_img_mm, Z_img_mm, image_fg, 'AlphaData', arg.transparencyMap);
        hold off;
    end

end

axis image

% color map handling
colormap(gca,my_map)
cbn = round(linspace(fgm,fgM,4),5);

% colorbar
if arg.showColorbar
    cb = colorbar('Ticks',[0, 0.33, 0.66, 1],'Ticklabels',...
        {num2str(cbn(1)), num2str(cbn(2)), num2str(cbn(3)), num2str(cbn(4))});
    set(cb,'YLim',[0 1])
    cb.Label.String = arg.ColorBarTitle;
end

% scale bar
if arg.scalebar>0    
    bar_pos = floor(0.9*max(Z_img_mm));
    hold on; line([1 2],[bar_pos bar_pos],'Color','white','LineWidth',5)
    text(1, bar_pos, '1 mm', 'Color', 'white',...
        'FontName', 'Helvetica', 'FontSize', 14,'verticalalignment','bottom')
end

% type set
ylabel(gca,'mm');
xlabel(gca,'mm');
ax = gca; ax.FontSize = arg.FontSize;

end