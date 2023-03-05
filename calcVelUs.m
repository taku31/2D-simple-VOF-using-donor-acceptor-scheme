function[UNs_new] = calcVelUs(UNs,NF, Pr, Uold, DX, Dtv, nx, ny, rho)

% 速度圧力補正：速度は液体要素のみ求める，粘性項，対流項，表面張力項を無視する

UNs_new = UNs;
ite_max = 1000;
eps  = 0.0001;

% 液体 NF = 0, 右液の垂直左表面 NF = 1, 左液の垂直右表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny
            
            if NF(i, j) == 2% 左液の垂直右表面ならば
                
                UNs_new(i, j) = UNs_new(i, j - 1);% 表面の流速は内部の流速を適用する。
                
            elseif  NF(i, j + 1) == 1 % 右液の垂直左表面ならば
                
                UNs_new(i, j) = UNs_new(i, j + 1);% 表面の流速は内部の流速を適用する。
                
            elseif  NF(i, j) == 5 && NF(i, j + 1) >= 1% 上表面かつその右側の要素が内部液体ではないならば。
                
                UNs_new(i, j) = UNs_new(i + 1, j );% 下の要素の速度を適用する。
                
            elseif NF(i, j) >= 1 && NF(i, j + 1) == 5% 内部液体ではないかつ右側が上表面ならば。？？？
                
                UNs_new(i, j) = UNs_new(i + 1, j);% 下の要素の速度を適用する。
                
            else
                
                UNs_new(i, j) = Uold(i, j) - (Pr(i, j + 1) - Pr(i, j)) / (DX * rho)  * Dtv;
                
            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(UNs_new - UNs)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    UNs = UNs_new;
    
end

end
