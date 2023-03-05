function[DTp_flag] = vis_contour(filename, timestep, time, DTp_flag, u, maxrange, minrange, fignum)

DTp_cnt = sum(DTp_flag(:, 2)) + 1;% “’B‚µ‚½ŠÔ‚ğŠÇ—‚·‚é•Ï”

if DTp_flag(DTp_cnt, 1) <= time
    
    f = figure(fignum);
    imagesc(u)
    colorbar
    caxis([minrange maxrange])
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(['VOF\_time = ', num2str(time, '%.3f [s]')]);
    axis equal; axis tight; axis on;
    
    % gif‚Å•Û‘¶
    frame = getframe(fignum);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if timestep == 1
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.0001, 'loopcount', inf);
    elseif rem(timestep, 1) == 0
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.0001, 'WriteMode', 'append');
    end
    
    % •Û‘¶‚µ‚½‚Éƒtƒ‰ƒO‚ğ—§‚Ä‚é
    DTp_flag(DTp_cnt, 2) = 1;
    
end

end