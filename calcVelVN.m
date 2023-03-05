function[VN_new] = calcVelVN(VN, NF, Pc, VNs, DX, nx, ny, rho,Dtv)

VN_new = VN;
ite_max = 1000;
eps  = 0.0001;

for k = 1 : ite_max
    
    for i = 2 : nx
        for j = 2 : ny + 1
            
            if   NF(i + 2 , j) >= 1 && NF(i + 1, j) == 8
                VN_new(i, j) = 0;
                
            elseif   NF(i + 1 , j) == 5
                
                VN_new(i, j) = 0;
                
            elseif   NF(i, j) == 2 && NF(i + 1 , j) >= 1
                
                VN_new(i, j) = VN_new(i, j - 1);
                
            elseif   NF(i + 1, j) == 2 && NF(i , j) >= 1
                
                VN_new(i, j) = VN_new(i, j - 1) ;
                
            elseif   NF(i, j) == 1 && NF(i + 1, j) >= 1
                
                VN_new(i, j) = VN_new(i, j + 1) ;
                
            elseif  NF(i + 1, j) == 1 && NF(i, j) >= 1
                
                VN_new(i, j) = VN_new(i, j + 1) ;
                
            else
                
                VN_new(i, j) = VNs(i, j) + (Pc(i, j) - Pc(i + 1, j)) / (DX * rho) * Dtv ;
                
            end
            
        end
    end
    
    % �c���̌v�Z
    error(k, 1) = max(max(abs(VN_new - VN)));
    if k ~= 1
        error(k, 1) = min(error(k, 1), error(k - 1, 1));
    end
    
    % ��������
    if  error(k, 1) < eps
        break;
    end
    
    % �l�̍X�V
    VN = VN_new;
    
end

end
