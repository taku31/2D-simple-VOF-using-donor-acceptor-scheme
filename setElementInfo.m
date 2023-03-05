function[NF_new] = setElementInfo(FF, NF, nx, ny)

NF_new = NF;
ite_max = 1000;
error = zeros(ite_max, 1);
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            if   FF(i, j) < 0.0001%VOF�l��0.001�ȉ��Ȃ�΋C�̂Ƃ݂Ȃ�
               
                NF_new(i, j) = 8;%�C�̔���
                
            else
                NF_new(i, j) = 0;%�t�̔���
            end
        end
    end
        
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if   NF_new(i, j) == 0%�t�̂Ȃ��
                
                if NF_new(i - 1, j) == 8% i,j ����̗v�f���A�C�̂Ȃ�Ώ�\��(5)�ł���B
                    
                    NF_new(i, j) = 5;%��\�ʔ���
                    
                elseif NF_new(i, j + 1 ) == 8% i,j �E�אڗv�f���C��(8)�Ȃ�΁A���t�̐����E�\��(5)�ł���B
                    
                    NF_new(i, j) = 2;%���t�̐����E�\��(5)
                    
                elseif NF_new(i, j - 1 ) == 8% i,j ���אڗv�f���C��(8)�Ȃ�΁A�E�t�̐������\��(1)�ł���B
                    
                    NF_new(i, j) = 1;%�E�t�̐������\��(1)
                    
                end
            end
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(NF_new - NF)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    NF = NF_new;
    
    
end

end