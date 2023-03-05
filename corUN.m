function[UNcor_new] = corUN(UNcor, NF3, UN2, nx, ny)

% UN ���x�␳ 4
% ���͂��C�̂Ȃ� u = 0�CNF = 1 �܂��� NF = 2 �Ȃ�Ύ��ʕۑ��̌v�Z�CNF = 5 �Ȃ�Ή��Ɠ��� u
% UNcor
% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9


UNcor_new = UNcor;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF3(i, j) == 8 && NF3(i - 1, j) >=8 && NF3(i + 1, j) >=8 && NF3(i, j - 1) >=8 && NF3(i, j + 1) >=8
                % ���̗v�f���C�̂ŏ����̎��͗v�f���C�̂��ǂ̎��͗����͂O�ɂ���B
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j) == 2
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j + 1) == 1
                
                UNcor_new(i, j) = 0;
                
            elseif  NF3(i, j) == 5 && NF3(i, j + 1) >= 1
                
                UNcor_new(i, j) = UNcor_new(i + 1, j);
                
            elseif NF3(i, j + 1) == 5 && NF3(i, j) >= 1
                
                UNcor_new(i, j) = UNcor_new(i + 1, j);
                
            else
                
                UNcor_new(i, j) = UN2(i, j);
                
            end
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(UNcor_new - UNcor)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    UNcor = UNcor_new;
    
end

end