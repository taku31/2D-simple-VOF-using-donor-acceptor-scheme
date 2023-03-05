function[UNcor_new] = corUN(UNcor, NF3, UN2, nx, ny)

% UN 速度補正 4
% 周囲が気体なら u = 0，NF = 1 または NF = 2 ならば質量保存の計算，NF = 5 ならば下と同じ u
% UNcor
% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9


UNcor_new = UNcor;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF3(i, j) == 8 && NF3(i - 1, j) >=8 && NF3(i + 1, j) >=8 && NF3(i, j - 1) >=8 && NF3(i, j + 1) >=8
                % その要素が気体で勝つその周囲要素が気体か壁の時は流速は０にする。
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j) == 2
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j + 1) == 1
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j) == 5 && NF3(i, j + 1) >= 1
                
                UNcor_new(i, j) = UNcor_new(i + 1, j);
                
            elseif NF3(i, j + 1) == 5 && NF3(i, j) >= 1
                
                UNcor_new(i, j) = UNcor_new(i + 1, j);
                
            else
                
                UNcor_new(i, j) = UN2(i, j);
                
            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(UNcor_new - UNcor)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    UNcor = UNcor_new;
    
end

end