function[VNs_new] = calcVelVs(VNs, NF, Pr, Vold, DX, Dtv, nx, ny, rho)


% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

VNs_new = VNs;
ite_max = 1000;
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx
        for j = 2 : ny + 1
            
            if  NF(i + 1, j) == 5% ���̗v�f����\�ʂȂ�΁H
                
                VNs_new(i, j) = VNs_new(i + 1, j);% ���̉��̗�����K�p����B
                
            elseif NF(i, j) == 2 && NF(i + 1, j) >= 1
                
                VNs_new(i, j) = VNs_new(i, j - 1);% ���t�̐����E�\�� ���A���̗v�f�������t�̈ȊO�Ȃ�΂��̍��̗�����K�p����B�Z
                
            elseif NF(i + 1, j) == 2 && NF(i , j) >= 1% ���̗v�f�����t�̐����E�\�� ���A���̗v�f�������t�̈ȊO�Ȃ�΁Z
                
                VNs_new(i, j) = VNs_new(i, j - 1);%���̍��̗�����K�p����B
                
            elseif NF(i, j) == 1 && NF(i + 1 , j) >= 1% ���̗v�f���E�t�̐������\�� ���A���̉��̗v�f�������t�̈ȊO�Ȃ��
                
                VNs_new(i, j) = VNs_new(i, j + 1);%���̉E�̗�����K�p����B
                
            elseif  NF(i + 1, j) == 1 && NF(i, j) >= 1% ���̉��̗v�f���E�t�̐������\�� ���A���̗v�f�������t�̈ȊO�Ȃ��
                
                VNs_new(i, j) = VNs_new(i, j + 1);%���̉E�̗�����K�p����B
                
            elseif NF(i, j) == 8
                
                VNs_new(i, j) = 0;
                
            else
                
                VNs_new(i, j) = Vold(i, j)  - (Pr(i + 1, j) - Pr(i, j)) * Dtv / (DX * rho)  + 9.8 * Dtv ;
                
            end
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(VNs_new - VNs)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    VNs = VNs_new;
    
end

end
