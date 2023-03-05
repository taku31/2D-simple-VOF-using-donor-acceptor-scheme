function[Pr_new] = calcPressure(NF, Pr, nx, ny, DX,FF2, rho)

% �O���� FF �����ł̉��̈��� P/r�F�\�ʗv�f�� FF = 0.5 �Ɖ��肵�Ē��S�� P/r = 0 �Ƃ���C
% �t�v�f P/r �͏�� P/r �� (�Ð���/���x) �����Z����
% (�Ð���/���xr) = (FFgh)
% �t�� NF = 0, �E�t�̐����\�� NF = 1, ���t�̐����\�� NF = 2, ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9

Pr_new = Pr;
ite_max = 1000;
error = zeros(ite_max, 1);
eps  = 0.0001;

for k = 1 : ite_max
  
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF(i, j) == 0% �t�̂ł���Έ��͂��v�Z
                
                Pr_new(i, j) = Pr_new(i - 1, j) + rho * FF2(i, j) * 9.8 * DX;
                
            elseif   NF(i, j) <= 2 && NF(i, j + 1) >= 8 && NF(i, j - 1) >= 8% �����\�ʗv�f�ł����Ă��A���ׂ��C���E�t�ʂł�����͈̂��͂��v�Z�B
                
                Pr_new(i, j) = Pr_new(i - 1, j) + rho * FF2(i, j) * 9.8 * DX;
                
            else % ��\�ʁE��C�v�f�E�ǂ��炩���t�̂̐����\�ʗv�f�ɂ͈��͂͂O�i��C���j�Ƃ���B
                
                Pr_new(i, j) = 0;
                
            end            
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(Pr_new - Pr)));
    if k ~= 1   
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
        
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    Pr = Pr_new;
    
end
