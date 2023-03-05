% ������
clearvars;

% �O���[�o���ϐ��錾
global nx ny

% �v�f�p�����[�^�[
nx = 30;% �v�f��
ny = 30;
DX = 1; % �v�f�� [m]
rho = 1;% ���̖��x [kg/m]

% �ϐ��錾
NF = zeros(nx + 2, ny + 2);
UNs = zeros(nx + 2, ny + 2);
VNs = zeros(nx + 2, ny + 2);
UN = zeros(nx + 2, ny + 2);
VN = zeros(nx + 2, ny + 2);
UN2 = zeros(nx + 2, ny + 2);
VN2 = zeros(nx + 2, ny + 2);
UNcor = zeros(nx + 2, ny + 2);
VNcor = zeros(nx + 2, ny + 2);
Uold = zeros(nx + 2, ny + 2);
Vold = zeros(nx + 2, ny + 2);
Ds = zeros(nx + 2, ny + 2);
PN = zeros(nx + 2, ny + 2);
Pc = zeros(nx + 2, ny + 2);
Pr = zeros(nx + 2, ny + 2);
FF = zeros(nx + 2, ny + 2);
FFnew = zeros(nx + 2, ny + 2);
FFcor= zeros(nx + 2, ny + 2);
FFcor2= zeros(nx + 2, ny + 2);
DFU = zeros(nx + 2, ny + 2);
DFV = zeros(nx + 2, ny + 2);
DFVcor= zeros(nx + 2, ny + 2);
DFVcor2= zeros(nx + 2, ny + 2);
UVtm= zeros(nx + 2, ny + 2);
time = 0 ;

% ���ԃp�����[�^�[
itemax = 100000;% �ő唽����
tmax = 20;% �v�Z�I������
Dt0 = 0.1;% �������ԍ���
Dtv = Dt0;% UV�v�Z�Ɏg��Dt [s]
cfl = 0.2;% �N�[������
Dtp = 0.01;% �|�X�g�����̊Ԋu
viscontf = transpose ([0 : Dtp : tmax; zeros(1, size(0 : Dtp : tmax, 2)) ]);% �R���^�[�����̐i�s���Ǘ�����t���O
visvectf = transpose ([0 : Dtp : tmax; zeros(1, size(0 : Dtp : tmax, 2)) ]);% �x�N�g�������̐i�s���Ǘ�����t���O

% �t�̌`�� FF �̏���������ݒ�
FF(end + 1 - round(end * 8 / 10) : end - 1, 2 : round(end * 3 / 10) ) = 1;%�t�̈ʒu���w��B
FFtotal0 = sum(sum(FF));% �����t�̗ʂ��v�Z

% �Z���̑������̏����ݒ�
% �t�� NF = 0, �E�t�̍������\�� NF = 1, ���t�̉E�����\�� NF = 2,
% ��\�� NF = 5, �C�� NF = 8�C�ǁC���y�ѓV�� NF =  9
NF(1, :) = 9;
NF(nx + 2, :) = 9;
NF(1 : nx + 2, 1) = 9;
NF(1 : nx + 2, ny + 2) = 9;
NF2 = NF;
NF3 = NF;

% ���Ԕ��W
for ite = 1 : itemax
    
    % �Z���������̍X�V
    NF = setElementInfo(FF, NF, nx, ny);
    
    % �������͂̌v�Z
    Pr = calcPressure(NF, Pr, nx, ny, DX, FF, rho);
    
    % ���̑��x�̌v�Z
    UNs = calcVelUs(UNs, NF, Pr, Uold, DX, Dtv, nx, ny, rho);
    VNs = calcVelVs(VNs, NF, Pr, Vold, DX, Dtv, nx, ny, rho);
    
    % �A���̎��̌덷���v�Z
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            Ds(i, j) =  rho *((UNs(i, j) - UNs(i, j - 1))/DX...
                + (VNs(i, j) - VNs(i - 1, j))/DX)/Dtv ;
        end
    end
    
    % ���͂̃|�A�\��������������
    Pc = solvePoissonEq(Pc, UNs, VNs,NF, DX, nx, ny, rho, Ds);
    
    % ���x�␳
    UN = calcVelUN(UN, VN, NF, Pc, UNs, DX, nx, ny, rho, Dtv);
    VN = calcVelVN(VN, NF, Pc, VNs, DX, nx, ny, rho, Dtv);
    
    % ����UV �̌v�Z�Ɏg�� Dt ����VOF�l�̈ڗ��v�Z�Ɏg��Dt�ɕύX
    UVN = abs(cat(1, UN, VN));
    
    % �N�[���������l�������^�C���X�e�b�v���v�Z
    if max(max(UVN(:)), min(UVN(:))) < 0.01 % �������\�����������
        Dta = 100;% �K���ɑ傫�Ȓl��ݒ肵�Ă���
    else
        Dta = cfl * DX / (max(max(UVN(:)) ,min(UVN(:))));% �N�[���������l������DT
    end
    Dtf = min(Dt0, Dta); % �����ݒ��dt�ƃN�[���������l������dt�̏������ق���p����B
    
    % Dtf�ɍ��킹�ė����E���͂��C��
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            UN2(i, j) = Uold(i, j) + (UN(i, j) - Uold(i, j)) * Dtf / Dtv;
            VN2(i, j) = Vold(i, j) + (VN(i, j) - Vold(i, j)) * Dtf / Dtv;
            
            %���͕␳
            if NF(i, j) == 0
                PN(i, j) = Pr(i,j) + Pc(i,j) / Dtf;
            else
                PN(i, j) = 0;
            end
            
        end
    end
    
    % �ڗ���DFU, DFV���v�Z
    DFU = calcDFU(DFU, DFV, NF, UN, Dtf, FF, DX, nx, ny);
    DFV = calcDFV(DFU, DFV, NF, VN, Dtf, FF, DX, nx, ny);
          
    % ���߂��ڗ���DFU, DFV��p����FF�̒l���X�V
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            FFnew(i, j) = DFU(i, j - 1) - DFU(i, j) +  DFV(i - 1, j) - DFV(i, j)  + FF(i, j) ;
        end
    end
    
    % �Z���������̍X�V
    NF2 = setElementInfo(FFnew, NF2, nx, ny);
    
    % DFV�␳�i1��ځj
    % �t���̃{�C�h�ƋC�̒��̉t�H����\�� NF = 5 �ɏW�߂邽�߂�DFV��␳����DFVcor�����߂�
    DFVcor = corDFV1(DFVcor, NF2, DFV, FFnew, nx, ny);
    
    % �␳�����ڗ���DFVcor��p����FF�̒l���X�V
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            FFcor(i, j) = DFU(i, j - 1) -  DFU(i, j) +  DFVcor(i - 1, j) - DFVcor(i, j) + FF(i, j);
            
        end
    end
    
    % DFV�␳�i2��ځj
    DFVcor2 = corDFV2(DFVcor2, DFU, DFVcor, NF, FF, nx, ny);
    
    % �␳�����ڗ���DFVcor2��p����FF�̒l���X�V
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            FFcor2(i, j) = DFU(i, j - 1) - DFU(i, j) + DFVcor2(i - 1, j) - DFVcor2(i, j) + FF(i, j);
            
        end
    end
    
    % �Z���������̍X�V
    NF3 = setElementInfo(FFcor2, NF3, nx, ny);
    
    % NF3�Ɋ�Â����x�␳
    UNcor = corUN(UNcor, NF3, UN2, nx, ny);
    VNcor = corVN(VNcor, NF3, VN2, nx, ny);
    
    % ���x�␳ UV �̍��v�l
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF3(i,j) == 0% �t�̕����̘A���̎����l����B
                
                UVtm(i, j) = UN(i, j - 1) - UN(i, j) + VN(i - 1, j) - VN(i, j);
                
            else
                
                UVtm(i, j) = 0;
                
            end
        end
    end
    UVtotal = sum(sum(UVtm));% �A���̎��̌덷
    UVmax = max(max(UVtm));% �ő�
    
    % �t�̂̑��ʂ��v�Z
    FFtotal = sum(sum(FFcor2));
    FFtotalp =  (FFtotal)/  FFtotal0 * 100;
    FFtotalps(ite, 1 : 3) = [ite, time, FFtotalp];
    
    % ����
    viscontf = vis_contour('vof.gif', ite, time, viscontf, FFcor2, 1, 0, 1);
    visvectf = vis_vector('velvec.gif', ite, time, visvectf, UNcor, VNcor, 2);
    vis_plot('totalvof.gif', time, tmax, ite, FFtotalps, 3)
    
    % ���Ԃ�i�߂�
    time = time + Dtf
    if time > tmax
        disp('�v�Z�I������tmax�ɓ��B�������ߌv�Z�I�����܂��B')
        break;
    end
    
    % �ߋ��̔z��Ɉړ�
    Dtv = Dtf;
    FF = FFcor2;
    Uold = UNcor;
    Vold = VNcor;
    
end
