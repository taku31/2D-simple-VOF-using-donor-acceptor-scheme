function[DFV_new] = calcDFV(DFU, DFV, NF, VN, Dtf, FF, DX, nx, ny)

% FF 移動量DFU,DFVの計算
% 但し，FFの移動量DFU,DFVの総和は上流のセルのFFを超えず（上流のセルのFFが負にならない），
% 自分のセルの空間部分 (1-FF) を超えない（自分のセルのFFが1を超えない）
% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

DFV_new = DFV;
ite_max = 10000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if   VN(i,j)>0 && NF(i, j) == 8
                
                DFV_new(i, j) =  0;
                
            elseif    NF(i + 1, j) == 0 && (NF(i, j) == 0 || VN(i,j) <= 0)
                
                DFV_new(i, j) =  - DFU(i + 1, j - 1) + DFU(i + 1, j) +  DFV_new(i + 1, j) + (1 - FF(i + 1, j));
                
            elseif    NF(i + 1, j) == 0
                
                DFV_new(i, j) =  min(- DFU(i + 1, j - 1) + DFU(i + 1, j) + DFV_new(i + 1, j) + (1 - FF(i + 1, j))...
                    , max(FF(i, j) + DFU(i, j - 1) - DFU(i, j) - DFV_new(i - 1, j) , 0));
                
            elseif   VN(i, j) >0 && NF(i, j) == 5
                
                DFV_new(i, j) = min( Dtf * VN(i, j) / DX * 1, FF(i, j) + DFU(i, j - 1) - DFU(i, j) - DFV_new(i - 1, j));
                
            elseif   VN(i, j) >0
                
                DFV_new(i, j) = Dtf * VN(i, j) / DX * FF(i, j);
                
            elseif   NF(i + 1, j) == 8
                
                DFV_new(i, j) = 0;
                
            elseif   NF(i + 1, j) == 5
                
                DFV_new(i, j) = - max(0, Dtf * abs(VN(i,j)) / DX * 1 - (1 - (FF(i + 1, j) + DFU(i + 1, j - 1) - DFU(i + 1, j)- DFV_new(i + 1, j) )));
                
            elseif   NF(i + 1, j) == 5
                
                DFV_new(i, j) = Dtf * VN(i,j) / DX * FF(i + 1, j) ;
                
            end
            
        end
        
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(DFV_new - DFV)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    DFV = DFV_new;
    
end

end
