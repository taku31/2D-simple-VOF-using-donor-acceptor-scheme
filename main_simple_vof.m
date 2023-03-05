% 初期化
clearvars;

% グローバル変数宣言
global nx ny

% 要素パラメーター
nx = 30;% 要素数
ny = 30;
DX = 1; % 要素幅 [m]
rho = 1;% 水の密度 [kg/m]

% 変数宣言
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

% 時間パラメーター
itemax = 100000;% 最大反復数
tmax = 20;% 計算終了時刻
Dt0 = 0.1;% 初期時間刻み
Dtv = Dt0;% UV計算に使うDt [s]
cfl = 0.2;% クーラン数
Dtp = 0.01;% ポスト処理の間隔
viscontf = transpose ([0 : Dtp : tmax; zeros(1, size(0 : Dtp : tmax, 2)) ]);% コンター可視化の進行を管理するフラグ
visvectf = transpose ([0 : Dtp : tmax; zeros(1, size(0 : Dtp : tmax, 2)) ]);% ベクトル可視化の進行を管理するフラグ

% 液体形状 FF の初期条件を設定
FF(end + 1 - round(end * 8 / 10) : end - 1, 2 : round(end * 3 / 10) ) = 1;%液体位置を指定。
FFtotal0 = sum(sum(FF));% 初期液体量を計算

% セルの属性情報の初期設定
% 液体 NF = 0, 右液の左垂直表面 NF = 1, 左液の右垂直表面 NF = 2,
% 上表面 NF = 5, 気体 NF = 8，壁，床及び天井 NF =  9
NF(1, :) = 9;
NF(nx + 2, :) = 9;
NF(1 : nx + 2, 1) = 9;
NF(1 : nx + 2, ny + 2) = 9;
NF2 = NF;
NF3 = NF;

% 時間発展
for ite = 1 : itemax
    
    % セル属性情報の更新
    NF = setElementInfo(FF, NF, nx, ny);
    
    % 初期圧力の計算
    Pr = calcPressure(NF, Pr, nx, ny, DX, FF, rho);
    
    % 仮の速度の計算
    UNs = calcVelUs(UNs, NF, Pr, Uold, DX, Dtv, nx, ny, rho);
    VNs = calcVelVs(VNs, NF, Pr, Vold, DX, Dtv, nx, ny, rho);
    
    % 連続の式の誤差を計算
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            Ds(i, j) =  rho *((UNs(i, j) - UNs(i, j - 1))/DX...
                + (VNs(i, j) - VNs(i - 1, j))/DX)/Dtv ;
        end
    end
    
    % 圧力のポアソン方程式を解く
    Pc = solvePoissonEq(Pc, UNs, VNs,NF, DX, nx, ny, rho, Ds);
    
    % 速度補正
    UN = calcVelUN(UN, VN, NF, Pc, UNs, DX, nx, ny, rho, Dtv);
    VN = calcVelVN(VN, NF, Pc, VNs, DX, nx, ny, rho, Dtv);
    
    % 流速UV の計算に使う Dt からVOF値の移流計算に使うDtに変更
    UVN = abs(cat(1, UN, VN));
    
    % クーラン数を考慮したタイムステップを計算
    if max(max(UVN(:)), min(UVN(:))) < 0.01 % 流速が十分小さければ
        Dta = 100;% 適当に大きな値を設定しておく
    else
        Dta = cfl * DX / (max(max(UVN(:)) ,min(UVN(:))));% クーラン数を考慮したDT
    end
    Dtf = min(Dt0, Dta); % 初期設定のdtとクーラン数を考慮したdtの小さいほうを用いる。
    
    % Dtfに合わせて流速・圧力を修正
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            UN2(i, j) = Uold(i, j) + (UN(i, j) - Uold(i, j)) * Dtf / Dtv;
            VN2(i, j) = Vold(i, j) + (VN(i, j) - Vold(i, j)) * Dtf / Dtv;
            
            %圧力補正
            if NF(i, j) == 0
                PN(i, j) = Pr(i,j) + Pc(i,j) / Dtf;
            else
                PN(i, j) = 0;
            end
            
        end
    end
    
    % 移流量DFU, DFVを計算
    DFU = calcDFU(DFU, DFV, NF, UN, Dtf, FF, DX, nx, ny);
    DFV = calcDFV(DFU, DFV, NF, VN, Dtf, FF, DX, nx, ny);
          
    % 求めた移流量DFU, DFVを用いてFFの値を更新
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            FFnew(i, j) = DFU(i, j - 1) - DFU(i, j) +  DFV(i - 1, j) - DFV(i, j)  + FF(i, j) ;
        end
    end
    
    % セル属性情報の更新
    NF2 = setElementInfo(FFnew, NF2, nx, ny);
    
    % DFV補正（1回目）
    % 液中のボイドと気体中の液滴を上表面 NF = 5 に集めるためにDFVを補正したDFVcorを求める
    DFVcor = corDFV1(DFVcor, NF2, DFV, FFnew, nx, ny);
    
    % 補正した移流量DFVcorを用いてFFの値を更新
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            FFcor(i, j) = DFU(i, j - 1) -  DFU(i, j) +  DFVcor(i - 1, j) - DFVcor(i, j) + FF(i, j);
            
        end
    end
    
    % DFV補正（2回目）
    DFVcor2 = corDFV2(DFVcor2, DFU, DFVcor, NF, FF, nx, ny);
    
    % 補正した移流量DFVcor2を用いてFFの値を更新
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            FFcor2(i, j) = DFU(i, j - 1) - DFU(i, j) + DFVcor2(i - 1, j) - DFVcor2(i, j) + FF(i, j);
            
        end
    end
    
    % セル属性情報の更新
    NF3 = setElementInfo(FFcor2, NF3, nx, ny);
    
    % NF3に基づき速度補正
    UNcor = corUN(UNcor, NF3, UN2, nx, ny);
    VNcor = corVN(VNcor, NF3, VN2, nx, ny);
    
    % 速度補正 UV の合計値
    for i = 2 : nx + 1
        for j = 2 : ny + 1
            
            if  NF3(i,j) == 0% 液体部分の連続の式を考える。
                
                UVtm(i, j) = UN(i, j - 1) - UN(i, j) + VN(i - 1, j) - VN(i, j);
                
            else
                
                UVtm(i, j) = 0;
                
            end
        end
    end
    UVtotal = sum(sum(UVtm));% 連続の式の誤差
    UVmax = max(max(UVtm));% 最大
    
    % 液体の総量を計算
    FFtotal = sum(sum(FFcor2));
    FFtotalp =  (FFtotal)/  FFtotal0 * 100;
    FFtotalps(ite, 1 : 3) = [ite, time, FFtotalp];
    
    % 可視化
    viscontf = vis_contour('vof.gif', ite, time, viscontf, FFcor2, 1, 0, 1);
    visvectf = vis_vector('velvec.gif', ite, time, visvectf, UNcor, VNcor, 2);
    vis_plot('totalvof.gif', time, tmax, ite, FFtotalps, 3)
    
    % 時間を進める
    time = time + Dtf
    if time > tmax
        disp('計算終了時間tmaxに到達したため計算終了します。')
        break;
    end
    
    % 過去の配列に移動
    Dtv = Dtf;
    FF = FFcor2;
    Uold = UNcor;
    Vold = VNcor;
    
end
