function[VNcor_new] = corVN(VNcor, NF, VN2, nx, ny)

VNcor_new = VNcor;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF(i, j) == 8 && NF(i - 1, j) >=8 && NF(i + 1, j - 1) >=8 && NF(i, j - 1) >=8 && NF(i, j + 1) >=8
                
                VNcor_new(i, j) = 0;
                                
            elseif  NF(i + 1, j) == 5
                
                VNcor_new(i, j) = 0;
                
            elseif NF(i, j) == 2 && NF(i + 1, j) >= 1
                
                VNcor_new(i, j) = VNcor_new(i, j - 1);

            elseif  NF(i + 1, j) == 2 && NF(i, j) >= 1
                
                VNcor_new(i, j) = VNcor_new(i, j - 1);
                
            elseif  NF(i, j + 1) == 1 && NF(i + 1, j + 1) >= 1
                
                VNcor_new(i, j) = VNcor_new(i, j + 1);
                
            elseif  NF(i + 1, j + 1) == 1 && NF(i, j + 1) >= 1
                
                VNcor_new(i, j) = VNcor_new(i, j + 1);
                
            else
                
                VNcor_new(i, j) = VN2(i, j);

            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(VNcor_new - VNcor)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    VNcor = VNcor_new;
    
end

end