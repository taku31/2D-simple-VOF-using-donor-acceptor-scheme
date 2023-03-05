function[DFV_new] = calcDFV(DFU, DFV, NF, VN, Dtf, FF, DX, nx, ny)

% FF �ړ���DFU,DFV�̌v�Z
% �A���CFF�̈ړ���DFU,DFV�̑��a�͏㗬�̃Z����FF�𒴂����i�㗬�̃Z����FF�����ɂȂ�Ȃ��j�C
% �����̃Z���̋�ԕ��� (1-FF) �𒴂��Ȃ��i�����̃Z����FF��1�𒴂��Ȃ��j
% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

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
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(DFV_new - DFV)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    DFV = DFV_new;
    
end

end
