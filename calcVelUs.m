function[UNs_new] = calcVelUs(UNs,NF, Pr, Uold, DX, Dtv, nx, ny, rho)

% ���x���͕␳�F���x�͉t�̗v�f�̂݋��߂�C�S�����C�Η����C�\�ʒ��͍��𖳎�����

UNs_new = UNs;
ite_max = 1000;
eps  = 0.0001;

% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny
            
            if NF(i, j) == 2% ���t�̐����E�\�ʂȂ��
                
                UNs_new(i, j) = UNs_new(i, j - 1);% �\�ʂ̗����͓����̗�����K�p����B
                
            elseif  NF(i, j + 1) == 1 % �E�t�̐������\�ʂȂ��
                
                UNs_new(i, j) = UNs_new(i, j + 1);% �\�ʂ̗����͓����̗�����K�p����B
                
            elseif  NF(i, j) == 5 && NF(i, j + 1) >= 1% ��\�ʂ����̉E���̗v�f�������t�̂ł͂Ȃ��Ȃ�΁B
                
                UNs_new(i, j) = UNs_new(i + 1, j );% ���̗v�f�̑��x��K�p����B
                
            elseif NF(i, j) >= 1 && NF(i, j + 1) == 5% �����t�̂ł͂Ȃ����E������\�ʂȂ�΁B�H�H�H
                
                UNs_new(i, j) = UNs_new(i + 1, j);% ���̗v�f�̑��x��K�p����B
                
            else
                
                UNs_new(i, j) = Uold(i, j) - (Pr(i, j + 1) - Pr(i, j)) / (DX * rho)  * Dtv;
                
            end
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(UNs_new - UNs)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    UNs = UNs_new;
    
end

end
