function[DFVcor_new] = corDFV1(DFVcor, NF, DFV, FFnew, nx, ny)

% DFV 補正 1
% 液中のボイドと気体中の液滴を上表面 NF = 5 に集めるためにDFVを補正したDFVcorを求めに行く。
% 液体 NF = 0, 右液の垂直表面 NF = 1, 左液の垂直表面 NF = 2, 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

DFVcor_new = DFVcor;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF(i, j) ==8% 気体ならば
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i - 1, j) - DFV(i - 1, j) + FFnew(i, j) ;
                
            elseif  NF(i + 1, j) == 0% 液体ならば
               
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i + 1, j) - DFV(i + 1, j) +  1 - FFnew(i + 1, j) ;
                
            elseif (NF(i + 1, j) == 2 || NF(i + 1, j) == 1)% 上の要素が 右液の垂直表面 NF = 1, 左液の垂直表面 NF = 2,ならば
               
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i + 1, j) - DFV(i + 1, j) ;
            
            else % そのNF2(i, j)は液体、左液の垂直表面、上表面のいずれか、NF2(i+1, j)は床および天井
                
                DFVcor_new(i, j) = DFV(i, j)  ;
                
            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(DFVcor_new - DFVcor)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    DFVcor = DFVcor_new;
    
end

end
