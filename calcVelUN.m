function[UN_new] = calcVelUN(UN, VN,NF, Pc,UNs,DX, nx, ny, rho, Dtv)

% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

UN_new = UN;
ite_max = 1000;
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny
            
            if   NF(i, j) == 2 % 左液の垂直右表面 NF = 2
                
                UN_new(i, j) = 0;
                
            elseif  NF(i, j + 1) == 1 %右どなり要素が 右液の垂直左表面 NF = 1ならば
                
                UN_new(i, j) = UN_new(i, j + 1) + VN(i, j + 1) - VN(i - 1, j + 1);
                
            elseif  NF(i, j) == 5 && NF(i, j + 1) >= 1
                
                UN_new(i, j) = UN_new(i + 1, j);
                
            elseif  NF(i, j + 1) == 5 && NF(i, j) >= 1
                
                UN_new(i, j) = UN_new(i + 1, j);
                
            else
                
                UN_new(i, j) = UNs(i, j)+(Pc(i, j) - Pc(i, j + 1)) / (DX*rho) * Dtv;
                
            end
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(UN_new - UN)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    UN = UN_new;
    
end

end
