function[NF_new] = setElementInfo(FF, NF, nx, ny)

NF_new = NF;
ite_max = 1000;
error = zeros(ite_max, 1);
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            if   FF(i, j) < 0.0001%VOF値が0.001以下ならば気体とみなす
               
                NF_new(i, j) = 8;%気体判定
                
            else
                NF_new(i, j) = 0;%液体判定
            end
        end
    end
        
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if   NF_new(i, j) == 0%液体ならば
                
                if NF_new(i - 1, j) == 8% i,j より上の要素が、気体ならば上表面(5)である。
                    
                    NF_new(i, j) = 5;%上表面判定
                    
                elseif NF_new(i, j + 1 ) == 8% i,j 右隣接要素が気体(8)ならば、左液の垂直右表面(5)である。
                    
                    NF_new(i, j) = 2;%左液の垂直右表面(5)
                    
                elseif NF_new(i, j - 1 ) == 8% i,j 左隣接要素が気体(8)ならば、右液の垂左直表面(1)である。
                    
                    NF_new(i, j) = 1;%右液の垂左直表面(1)
                    
                end
            end
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(NF_new - NF)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    NF = NF_new;
    
    
end

end