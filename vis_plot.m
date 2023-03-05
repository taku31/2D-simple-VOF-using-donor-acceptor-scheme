function[] = vis_plot(filename, time, tmax, timestep, FFtotals, fignum)

figure(fignum);
h = plot(FFtotals(:, 2), FFtotals(:, 3)) ;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
title(['TotalVOF\_time = ', num2str(time, '%.3f [s]')]);
axis tight; axis on;
xlim([0 tmax]);
ylim([0 100 * 1.1]);
xlabel('Time [s]')
ylabel('TotalVOF [%]')


% gif‚Å•Û‘¶
frame = getframe(fignum);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if timestep == 1
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'loopcount', inf);
elseif rem(timestep, 1) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
end

end