function[DFU_new] = calcDFU(DFU, DFV, NF, UN, Dtf, FF, DX, nx, ny)

% FF 移動量DFU,DFVの計算
% 但し，FFの移動量DFU,DFVの総和は上流のセルのFFを超えず（上流のセルのFFが負にならない），
% 自分のセルの空間部分 (1-FF) を超えない（自分のセルのFFが1を超えない）
% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

DFU_new = DFU;
ite_max = 10000;
eps  = 1e-6;

for k = 1 : ite_max
    
    % DFUの計算
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            % 流速が正の場合（F(i, j) がドナーセルになる場合）
            
            if   UN(i, j)>0 && NF(i, j) == 8 %その要素の流速が正で、気体ならば、
                
                DFU_new(i, j) = 0 ;
                
            elseif   UN(i, j)>0 && NF(i, j) == 1
                
                DFU_new(i, j) = min(Dtf * UN(i,j) / DX * 1, FF(i, j) +  DFV(i - 1, j)-  DFV(i, j) +  DFU_new(i, j - 1)) ;
                
            elseif   UN(i, j)>=0 && NF(i, j) == 2
                
                DFU_new(i, j) = max(0, (Dtf * UN(i,j) / DX * 1-(1 -(FF(i, j) +  DFV(i - 1, j)-  DFV(i, j) +  DFU_new(i, j - 1) ))));
                
            elseif   UN(i, j)>=0
                
                DFU_new(i, j) = Dtf * UN(i,j) / DX * FF(i, j);
                
                % 流速が負の場合（F(i, j) がアクセプターセルになる場合）
                
            elseif   NF(i, j + 1) == 8
                
                DFU_new(i, j) = 0;
                
            elseif   NF(i, j + 1) == 2
                
                DFU_new(i, j) = - min(Dtf * abs(UN(i,j)) / DX * 1, FF(i, j + 1) + DFV(i - 1, j + 1) - DFV(i, j + 1) - DFU_new(i, j + 1) );
                
            elseif   NF(i, j + 1) == 1
                
                DFU_new(i, j) = - max(0, Dtf * abs(UN(i,j)) / DX * 1 - (FF(i,j + 1) + DFV(i - 1, j + 1) - DFV(i, j + 1) - DFU_new(i, j + 1)) );
                
            else
                
                DFU_new(i, j) =  Dtf * UN(i,j) / DX * FF(i, j + 1);
                
            end
            
        end
    end
    
 
    % 残差の計算
    error(k, 1) = max(max(abs(DFU_new - DFU)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    DFU = DFU_new;
    
end

end
