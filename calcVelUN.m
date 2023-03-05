function[UN_new] = calcVelUN(UN, VN,NF, Pc,UNs,DX, nx, ny, rho, Dtv)

% �t�� NF = 0, �E�t�̐������\�� NF = 1, ���t�̐����E�\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

UN_new = UN;
ite_max = 1000;
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny
            
            if   NF(i, j) == 2 % ���t�̐����E�\�� NF = 2
                
                UN_new(i, j) = 0;
                
            elseif  NF(i, j + 1) == 1 %�E�ǂȂ�v�f�� �E�t�̐������\�� NF = 1�Ȃ��
                
                UN_new(i, j) = UN_new(i, j + 1) + VN(i, j + 1) - VN(i - 1, j + 1);
                
            elseif  NF(i, j) == 5 && NF(i, j + 1) >= 1
                
                UN_new(i, j) = UN_new(i + 1, j);
                
            elseif  NF(i, j + 1) == 5 && NF(i, j) >= 1
                
                UN_new(i, j) = UN_new(i + 1, j);
                
            else
                
                UN_new(i, j) = UNs(i, j)+(Pc(i, j) - Pc(i, j + 1)) / (DX*rho) * Dtv;
                
            end
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(UN_new - UN)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    UN = UN_new;
    
end

end
