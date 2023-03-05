function[DFVcor_new] = corDFV1(DFVcor, NF, DFV, FFnew, nx, ny)

% DFV �␳ 1
% �t���̃{�C�h�ƋC�̒��̉t�H����\�� NF = 5 �ɏW�߂邽�߂�DFV��␳����DFVcor�����߂ɍs���B
% �t�� NF = 0, �E�t�̐����\�� NF = 1, ���t�̐����\�� NF = 2, ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

DFVcor_new = DFVcor;
ite_max = 1000;
eps  = 1e-6;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF(i, j) ==8% �C�̂Ȃ��
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i - 1, j) - DFV(i - 1, j) + FFnew(i, j) ;
                
            elseif  NF(i + 1, j) == 0% �t�̂Ȃ��
               
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i + 1, j) - DFV(i + 1, j) +  1 - FFnew(i + 1, j) ;
                
            elseif (NF(i + 1, j) == 2 || NF(i + 1, j) == 1)% ��̗v�f�� �E�t�̐����\�� NF = 1, ���t�̐����\�� NF = 2,�Ȃ��
               
                DFVcor_new(i, j) = DFV(i, j) + DFVcor_new(i + 1, j) - DFV(i + 1, j) ;
            
            else % ����NF2(i, j)�͉t�́A���t�̐����\�ʁA��\�ʂ̂����ꂩ�ANF2(i+1, j)�͏�����ѓV��
                
                DFVcor_new(i, j) = DFV(i, j)  ;
                
            end
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(DFVcor_new - DFVcor)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    DFVcor = DFVcor_new;
    
end

end
