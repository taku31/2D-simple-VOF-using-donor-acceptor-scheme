function[DFU_new] = calcDFU(DFU, DFV, NF, UN, Dtf, FF, DX, nx, ny)

% FF �ړ���DFU,DFV�̌v�Z
% �A���CFF�̈ړ���DFU,DFV�̑��a�͏㗬�̃Z����FF�𒴂����i�㗬�̃Z����FF�����ɂȂ�Ȃ��j�C
% �����̃Z���̋�ԕ��� (1-FF) �𒴂��Ȃ��i�����̃Z����FF��1�𒴂��Ȃ��j
% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

DFU_new = DFU;
ite_max = 10000;
eps  = 1e-6;

for k = 1 : ite_max
    
    % DFU�̌v�Z
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            % ���������̏ꍇ�iF(i, j) ���h�i�[�Z���ɂȂ�ꍇ�j
            
            if   UN(i, j)>0 && NF(i, j) == 8 %���̗v�f�̗��������ŁA�C�̂Ȃ�΁A
                
                DFU_new(i, j) = 0 ;
                
            elseif   UN(i, j)>0 && NF(i, j) == 1
                
                DFU_new(i, j) = min(Dtf * UN(i,j) / DX * 1, FF(i, j) +  DFV(i - 1, j)-  DFV(i, j) +  DFU_new(i, j - 1)) ;
                
            elseif   UN(i, j)>=0 && NF(i, j) == 2
                
                DFU_new(i, j) = max(0, (Dtf * UN(i,j) / DX * 1-(1 -(FF(i, j) +  DFV(i - 1, j)-  DFV(i, j) +  DFU_new(i, j - 1) ))));
                
            elseif   UN(i, j)>=0
                
                DFU_new(i, j) = Dtf * UN(i,j) / DX * FF(i, j);
                
                % ���������̏ꍇ�iF(i, j) ���A�N�Z�v�^�[�Z���ɂȂ�ꍇ�j
                
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
    
 
    % �c���̌v�Z
    error(k, 1) = max(max(abs(DFU_new - DFU)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    DFU = DFU_new;
    
end

end
