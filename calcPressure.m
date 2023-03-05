function[Pr_new] = calcPressure(NF, Pr, nx, ny, DX,FF2, rho)

% 前時刻 FF 条件での仮の圧力 P/r：表面要素は FF = 0.5 と仮定して中心の P/r = 0 とする，
% 液要素 P/r は上の P/r に (静水圧/密度) を加算する
% (静水圧/密度r) = (FFgh)
% 液体 NF = 0, 右液の垂直表面 NF = 1, 左液の垂直表面 NF = 2, 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9

Pr_new = Pr;
ite_max = 1000;
error = zeros(ite_max, 1);
eps  = 0.0001;

for k = 1 : ite_max
  
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF(i, j) == 0% 液体であれば圧力を計算
                
                Pr_new(i, j) = Pr_new(i - 1, j) + rho * FF2(i, j) * 9.8 * DX;
                
            elseif   NF(i, j) <= 2 && NF(i, j + 1) >= 8 && NF(i, j - 1) >= 8% 垂直表面要素であっても、両隣が気相・液面であるものは圧力を計算。
                
                Pr_new(i, j) = Pr_new(i - 1, j) + rho * FF2(i, j) * 9.8 * DX;
                
            else % 上表面・空気要素・どちらかが液体の垂直表面要素には圧力は０（大気圧）とする。
                
                Pr_new(i, j) = 0;
                
            end            
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(Pr_new - Pr)));
    if k ~= 1   
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
        
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    Pr = Pr_new;
    
end
