function[DTp_flag] = vis_vector(filename, timestep, time, DTp_flag, u, v, fignum)

% グローバル変数呼び出し
global nx ny

DTp_cnt = sum(DTp_flag(:, 2)) + 1;% 到達した時間を管理する変数

if DTp_flag(DTp_cnt, 1) <= time
    
    figure(fignum);
    h = quiver(flipud(u),-flipud(v)) ;
    h.Color = 'r';
    h.AutoScaleFactor= 1.5;
    
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(['velvec\_time = ', num2str(time, '%.3f [s]')]);
    axis equal; axis tight; axis on;
    xlim([1 nx + 2]);
    ylim([1 ny + 2]);
    
    % gifで保存
    frame = getframe(fignum);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if timestep == 1
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'loopcount', inf);
    elseif rem(timestep, 1) == 0
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
    end
    
    % 保存した時刻にフラグを立てる
    DTp_flag(DTp_cnt, 2) = 1;
    
end

end