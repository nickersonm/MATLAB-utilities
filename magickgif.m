% Use ImageMagick to create an animated gif
function magickgif(imgarray, outFile, varargin)
    delay = 100/4;
    if ~isempty(varargin) > 0
        delay = 100/varargin{1};    % FPS input
    end
    
    mkdir("./tmpmagick");
    for i = 1:length(imgarray)
        img = imgarray(i);
        if isa(img, 'struct')   % Likely a frame
            img = frame2im(img);
        end
        if isa(img, 'cell')
            img = img{1};
        end
        
        imwrite(img, sprintf('./tmpmagick/%03i.png', i));
    end
    system( sprintf("magick convert -delay %i ./tmpmagick/*.png -loop 0 -coalesce -layers Optimize +map %s", delay, outFile));
    rmdir './tmpmagick' s;
end
