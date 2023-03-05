function[DFVcor2_new] = corDFV2(DFVcor2, DFU, DFVcor, NF, FF_old, nx, ny)

%DFV 補正 2
%上の FF < 0 の場合に下の DFV を上の FF = 0 となるように変える，
%上の FF が下の FF より大なら FF を平均値にする，
%但し天井直下の一番上の層は除く，液体内部で FF の多い分を吸収できなければ気体中に FF を移動させる

DFVcor2_new = DFVcor2;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx
        for j = 2 : ny + 1
            
            %上の FF < 0 の場合に下の DFV を上の FF = 0 となるように変える，
            if ( DFU(i, j - 1) - DFU(i, j) + DFVcor(i - 1, j)  - DFVcor(i, j)  + FF_old(i, j) < 0 ...
                    || DFU(i, j - 1) - DFU(i, j) + DFVcor2_new(i - 1, j) - DFVcor2_new(i, j) + FF_old(i, j) - 0.01 < 0)
                
                DFVcor2_new(i, j)= min(DFU(i, j - 1)     - DFU(i, j)     + DFVcor2_new(i - 1, j) + FF_old(i, j), ...
                    1 - (DFU(i + 1, j - 1) - DFU(i + 1, j) - DFVcor2_new(i + 1, j) + FF_old(i + 1, j)));
                
            elseif  NF(i - 1, j) == 9
                
                DFVcor2_new(i, j) = DFVcor(i, j);
                                
            elseif ((DFU(i, j - 1)     - DFU(i, j)     + DFVcor(i - 1, j)  - DFVcor(i, j)      + FF_old(i, j)) > ...
                    (DFU(i + 1, j - 1) - DFU(i + 1, j) + DFVcor(i, j)      - DFVcor(i + 1, j)  + FF_old(i + 1, j) + 0.01)...
                    || (DFU(i, j - 1)     - DFU(i, j)     + DFVcor2_new(i - 1, j) - DFVcor2_new(i, j)     + FF_old(i, j)) > ...
                    (DFU(i + 1, j - 1) - DFU(i + 1, j) + DFVcor2_new(i, j)     - DFVcor2_new(i + 1, j) + FF_old(i + 1, j) + 0.01))
                
                DFVcor2_new(i, j) = min(((DFU(i, j - 1)     - DFU(i, j)     + DFVcor2_new(i - 1, j) + FF_old(i, j))...
                    - (DFU(i + 1, j - 1) - DFU(i + 1, j) - DFVcor2_new(i + 1, j) + FF_old(i + 1, j)))/2, ...
                    1 - (DFU(i + 1, j - 1) - DFU(i + 1, j) - DFVcor2_new(i + 1, j) + FF_old(i + 1, j)));
                
            elseif DFU(i, j - 1)     - DFU(i, j)     + DFVcor2_new(i - 1, j) - DFVcor2_new(i, j)     + FF_old(i, j) < 1 ...
                    && DFU(i + 1, j - 1) - DFU(i + 1, j) + DFVcor2_new(i, j)     - DFVcor2_new(i + 1, j) + FF_old(i + 1, j) > 1
                
                DFVcor2_new(i, j) = 1 - (DFU(i + 1, j - 1) - DFU(i + 1, j)                     - DFVcor2_new(i + 1, j) + FF_old(i + 1, j));
               
            else
                
                DFVcor2_new(i, j) = DFVcor(i, j);
            end
            
        end
    end
    
    % 残差の計算
    error(k, 1) = max(max(abs(DFVcor2_new - DFVcor2)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % 収束判定
    if  error(k, 1) < eps
        break;
    end
    
    % 値の更新
    DFVcor2 = DFVcor2_new;   
    
end

end
