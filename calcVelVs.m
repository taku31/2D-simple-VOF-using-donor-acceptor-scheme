function[VNs_new] = calcVelVs(VNs, NF, Pr, Vold, DX, Dtv, nx, ny, rho)


% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

VNs_new = VNs;
ite_max = 1000;
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx
        for j = 2 : ny + 1
            
            if  NF(i + 1, j) == 5% 下の要素が上表面ならば？
                
                VNs_new(i, j) = VNs_new(i + 1, j);% その下の流速を適用する。
                
            elseif NF(i, j) == 2 && NF(i + 1, j) >= 1
                
                VNs_new(i, j) = VNs_new(i, j - 1);% 左液の垂直右表面 かつ、下の要素が内部液体以外ならばその左の流速を適用する。〇
                
            elseif NF(i + 1, j) == 2 && NF(i , j) >= 1% 下の要素が左液の垂直右表面 かつ、その要素が内部液体以外ならば〇
                
                VNs_new(i, j) = VNs_new(i, j - 1);%その左の流速を適用する。
                
            elseif NF(i, j) == 1 && NF(i + 1 , j) >= 1% その要素が右液の垂直左表面 かつ、その下の要素が内部液体以外ならば
                
                VNs_new(i, j) = VNs_new(i, j + 1);%その右の流速を適用する。
                
            elseif  NF(i + 1, j) == 1 && NF(i, j) >= 1% その下の要素が右液の垂直左表面 かつ、その要素が内部液体以外ならば
                
                VNs_new(i, j) = VNs_new(i, j + 1);%その右の流速を適用する。
                
            elseif NF(i, j) == 8
                
                VNs_new(i, j) = 0;
                
            else
                
                VNs_new(i, j) = Vold(i, j)  - (Pr(i + 1, j) - Pr(i, j)) * Dtv / (DX * rho)  + 9.8 * Dtv ;
                
            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(VNs_new - VNs)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    VNs = VNs_new;
    
end

end
