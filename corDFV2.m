function[DFVcor2_new] = corDFV2(DFVcor2, DFU, DFVcor, NF, FF_old, nx, ny)

%DFV �␳ 2
%��� FF < 0 �̏ꍇ�ɉ��� DFV ����� FF = 0 �ƂȂ�悤�ɕς���C
%��� FF ������ FF ����Ȃ� FF �𕽋ϒl�ɂ���C
%�A���V�䒼���̈�ԏ�̑w�͏����C�t�̓����� FF �̑��������z���ł��Ȃ���΋C�̒��� FF ���ړ�������

DFVcor2_new = DFVcor2;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx
        for j = 2 : ny + 1
            
            %��� FF < 0 �̏ꍇ�ɉ��� DFV ����� FF = 0 �ƂȂ�悤�ɕς���C
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
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(DFVcor2_new - DFVcor2)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    DFVcor2 = DFVcor2_new;   
    
end

end
