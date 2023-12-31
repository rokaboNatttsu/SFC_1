using StatsPlots

#   関数定義
function btw(X::Float64, Y::Float64, Z::Float64)
    if X > Z    #   不正入力に対してはエラーを吐く値を返す
        return nothing
    end
    return max(X, min(Y, Z))
end
function cal_estimated_value(X, Xe, t, T, λe)
    tm1 = (T + t - 1 - 1) % T + 1
    tm2 = (T + t - 2 - 1) % T + 1
    tm3 = (T + t - 3 - 1) % T + 1
    tm4 = (T + t - 4 - 1) % T + 1
    return λe*X[tm1] + (1 - λe)*Xe[tm1] + (X[tm1] - X[tm2])/2 + (X[tm2] - X[tm3])/3 + (X[tm3] - X[tm4])/6
end

#   外生変数
T = 200 + 1     #   > 100
m = 0.1 #   税も含めた総コストに対する価格マークアップ
ut, ugt = 0.75, 0.7    #   目標資本稼働率
Cg0, SS0 = 50.0, 10.0
ig, pe, i = 0.01, 1.0, 0.02

#   パラメータ
α1, α2, α3, α4, α5 = 0.9, 0.05, 0.1, 0.05, 0.0
βk, βc, βgk, βgc, βg1, βg2, βg3 = 0.4, 0.4, 0.4, 4.0, α5, 0.5, 0.3
β0, β2, β3, β4, β5, β6 = 0.1, 0.2, 0.1, 0.05, 1.0, 1.0
β1 = (βk*ut - β0)/(βk*ut)
γ1, γ2, γ3 = 0.02, 0.02, 0.02
γc, γk, γg = 1.0, 2.0, 2.0
δv, δiw, δii, δem, δei, δef, δf  = 0.1, 0.2, 0.2, 0.01, 0.02, 0.015, 0.2
ϵ, η, ζ = 0.65, 5.0, 0.01
λ00, λ10, λ20, λ01, λ11, λ21 = 0.4, 0.4, 0.2, -1.0, 1.0, 0.0
λ21 = -(λ01 + λ11)
λ20 = 1 - (λ00 + λ10)
λe, λp = 0.5, 0.75

#   変数定義
p, pk = zeros(T), zeros(T), zeros(T)
u, uc, uk, ug, ugk, ugc, ue = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
UCe, UC = zeros(T), zeros(T)
I, Ig, Igg, Igf, Ie, Igfe = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
K, Kg = zeros(T), zeros(T)
Cw, Ci, Cg, C, Ce = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
SS, C_demand = zeros(T), zeros(T)
ΔIN, IN, S, INe, Se = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
YDw, YDi, YDwe, YDie = zeros(T), zeros(T), zeros(T), zeros(T)
NWw, NWi, NWwe, NWie = zeros(T), zeros(T), zeros(T), zeros(T)
Tve, Tfe, Tefe = zeros(T), zeros(T), zeros(T)
Tv, Tf = zeros(T), zeros(T)
Tiw, Tii, Ti = zeros(T), zeros(T), zeros(T)
Tew, Tei, Tef, Te = zeros(T), zeros(T), zeros(T), zeros(T)
Nc, Nk, Ng = zeros(T), zeros(T), zeros(T)
Πi, Πf, Πb, Π = zeros(T), zeros(T), zeros(T), zeros(T)
Mw, Mi, Mf, M = zeros(T), zeros(T), zeros(T), zeros(T)
ΔMw, ΔMi, ΔMf, ΔM = zeros(T), zeros(T), zeros(T), zeros(T)
WN, WbNb, Wf, Wg = zeros(T), zeros(T), zeros(T), zeros(T)
GBi, GBb, GB = zeros(T), zeros(T), zeros(T)
ΔGBi, ΔGBb, ΔGB = zeros(T), zeros(T), zeros(T)
Lw, Li, Lf, L = zeros(T), zeros(T), zeros(T), zeros(T)
ΔLw, ΔLi, ΔLf, ΔL = zeros(T), zeros(T), zeros(T), zeros(T)
Ei, Eb, E = zeros(T), zeros(T), zeros(T)
ΔEi, ΔEb, ΔE = zeros(T), zeros(T), zeros(T)
GS = zeros(T)
NLw, NLi, NLf, NLg = zeros(T), zeros(T), zeros(T), zeros(T)
NWw, NWi, NWf, NWg = zeros(T), zeros(T), zeros(T), zeros(T)
H, ΔH = zeros(T), zeros(T)

#   初期値設定
u[end], ue[end] = ut, ut
Nc[end], Nk[end], Ng[end] = 1200.0, 1200.0, 180.0
UC[end], p[end], pk[end] = 0.2, 0.22, 0.22
I[end], Ie[end] = 1.0, 1.0
K[end], Kg[end] = 1800.0, 300.0
C[end], Ce[end], INe[end], IN[end], S[end] = 350.0, 350.0, 70.0, 70.0, 350.0
C_demand[end] = C[end]
YDwe[end], YDw[end], YDie[end], YDi[end] = 55.0, 55.0, 2.0, 2.0
NWwe[end], NWw[end], NWie[end], NWi[end] = 500.0, 500.0, 110.0, 110.0
Π[end] = 1.0
Wf[end], Wg[end], WN[end] = 0.01, 0.01, 6.0
Igf[end], Igfe[end], Πi[end], GBi[end], Tiw[end], Tii[end] = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
Ei[end] = 10.0
IN[1], Ei[1], E[1] = β2*C[end], Ei[end], Ei[end]
Mf[1], Mi[1], Mw[1] = 2500.0, 100.0, 500.0
M[1], H[1] = sum([Mf[1], Mi[1], Mw[1]]), sum([Mf[1], Mi[1], Mw[1]])
K[1], Kg[1] = K[end], Kg[end]
NWw[1] = Mw[1]
NWi[1] = pe*Ei[1] + Mi[1]
NWf[1] = pk[end]*K[1] + p[end]*IN[1] - pe*E[1] + Mf[1]
NWg[1] = p[end]*Kg[1] - H[1]
Wf[1] = Wf[end]
    #   ストックの変数の初期値は必ず会計的一貫性を持つように設定すること。未確認

#   定常状態に至るまでシミュレーション
for t = 1:T-1
    tm1 = (T + t - 1 - 1) % T + 1
    tm2 = (T + t - 2 - 1) % T + 1
    tm3 = (T + t - 3 - 1) % T + 1
    tm4 = (T + t - 4 - 1) % T + 1
    #   期待値計算
    ue[t] = cal_estimated_value(u, ue, t, T, λe)
    Ie[t] = cal_estimated_value(I, Ie, t, T, λe)
    Ce[t] = cal_estimated_value(C, Ce, t, T, λe)
    Se[t] = cal_estimated_value(S, Se, t, T, λe)
    INe[t] = cal_estimated_value(IN, INe, t, T, λe)
    Igfe[t] = cal_estimated_value(Igf, Igfe, t, T, λe)
    YDwe[t] = cal_estimated_value(YDw, YDwe, t, T, λe)
    YDie[t] = cal_estimated_value(YDi, YDie, t, T, λe)
    Tve[t] = cal_estimated_value(Tv, Tve, t, T, λe)
    Tfe[t] = cal_estimated_value(Tf, Tfe, t, T, λe)
    Tefe[t] = cal_estimated_value(Tef, Tefe, t, T, λe)
    NWwe[t] = cal_estimated_value(NWw, NWwe, t, T, λe)
    NWie[t] = cal_estimated_value(NWi, NWie, t, T, λe)
    #   雇用者数設定
    Nc[t] = btw((1-β4)*Nc[tm1], γc*β1*K[t], (1+β3)*Nc[tm1])    #   行動方程式
    Nk[t] = btw((1-β4)*Nk[tm1], γk*(1-β1)*K[t], (1+β3)*Nk[tm1])    #   行動方程式
    Ng[t] = btw((1-β4)*Ng[tm1], γg*βg3*Kg[t], (1+β3)*Ng[tm1])  #   行動方程式
    #   価格設定
    Wf[t] = Wf[tm1]*((1-γ1) + γ1*ue[t]/ut)  #   行動方程式
    Wg[t] = Wg[tm1]*((1-γ3) + γ3*Wf[tm1]/Wg[tm1])   #   行動方程式
    UCe[t] = (UC[tm1]*IN[t] + (Wf[t]*(Nc[t] + Nk[t]) + Tve[t] + Tfe[t] + Tefe[t]))/(Se[t] + Igfe[t] + IN[t]) #   行動方程式
    p[t] = λp*p[tm1] + (1 - λp)*(1.0 + m)*UCe[t]  #   行動方程式
    pk[t] = γk*Wf[t]/(βk*(1.0 - β1))   #   行動方程式 今から同じ量の資本を再取得するために必要な賃金の合計、といったニュアンス
    #   生産活動の意思決定
    Igf[t] = min((β0 + βg1)*Kg[t]/(1.0 + βg2), βc*Nc[t]/γc, βc*β1*K[t])    #   行動方程式
    Igg[t] = min((βg2*βg1 - β0)*Kg[t]/(1.0 + βg2), βgk*Ng[t]/γg, βgk*βg3*Kg[t] - β0*Kg[t])  #   行動方程式
    Ig[t] = Igf[t] + Igg[t] #   水平一貫性の会計恒等式
    S[t] = btw(0.0, Ce[t] - IN[t] + β2*Ce[t], max(0.0, min(βc*Nc[t]/γc - Igfe[t], βc*β1*K[t] - Igfe[t], βgc*(1.0 - βg3)*Kg[t])))   #   行動方程式
    I[t] = -β0*K[t] + btw(0.0, ((Se[t] + Igfe[t])/(βc*β1*K[t]) - ut)*K[t] + β0*K[t], min(βk*(1.0 - β1)*K[t], βk*Nk[t]/γk))#, max(0.0, ue[t]/ut*(β5*Mf[t] + β6*NLf[tm1]))))  #   行動方程式
    uc[t] = (S[t] + Igf[t])/(βc*β1*K[t])  #   定義式
    uk[t] = (I[t] + β0*K[t])/(βk*(1-β1)*K[t])   #   定義式
    u[t] = (S[t] + Igf[t])/(βc*K[t]) + (I[t] + β0*K[t])/(βk*K[t]) #   定義式
    ugk[t] = (Igg[t] + β0*Kg[t])/(βgk*βg3*Kg[t])    #   定義式
    ugc[t] = S[t]/(βgc*(1.0 - βg3)*Kg[t]) #   定義式
    ug[t] = (1 - βg3)*ugc[t] + βg3*ugk[t]   #   定義式
    K[t+1] = K[t] + I[t]    #   ストックとフローの整合性の会計恒等式
    Kg[t+1] = Kg[t] + Ig[t]    #   ストックとフローの整合性の会計恒等式
    #   支出の意思決定
    Cw[t] = (α1*YDwe[t] + α2*NWwe[t])/p[t]  #   行動方程式
    Ci[t] = (α3*YDie[t] + α4*NWie[t])/p[t]  #   行動方程式
    Cg[t] = Cg0*(1.0 + α5)^t  #   行動方程式
    C_demand[t] = Cw[t] + Ci[t] + Cg[t]
    C[t] = min(IN[t] + S[t], C_demand[t])    #   水平一貫性の会計恒等式
    if C_demand[t] > C[t]
        tmp = C[t]/C_demand[t]
        Cw[t] *= tmp
        Ci[t] *= tmp
        Cg[t] *= tmp
        println(t, " C_demand > C")
    end
    SS[t] = SS0*Kg[t]/Kg[1]  #   行動方程式
    IN[t+1] = IN[t] + S[t] - C[t]   #   行動方程式
    #   税額の決定  （循環参照を避けるため、仕方なく、１期前の値から計算する。税額は後から確定するんだから問題ないってな言い訳もできるっちゃできる）
    Tv[t] = max(0.0, δv*(p[tm1]*(C[tm1] + ΔIN[tm1] + Igf[tm1]) + pk[tm1]*I[tm1]) - i*Lf[t])   #   行動方程式
    Tiw[t] = δiw*WN[tm1]  #   行動方程式
    Tii[t] = δii*(Πi[tm1] + ig*GBi[tm1])    #   行動方程式
    Ti[t] = Tiw[t] + Tii[t] #   水平一貫性の会計恒等式
    Tew[t] = δem*max(0.0, Mw[tm1] - Lw[tm1])  #   行動方程式
    Tei[t] = δem*max(0.0, Mi[tm1] - Li[tm1]) + δei*(pe*Ei[tm1] + GBi[tm1]) #   行動方程式
    Tef[t] = δef*K[tm1] #   行動方程式
    Te[t] = Tew[t] + Tei[t] + Tef[t]    #   水平一貫性の会計恒等式
    Tf[t] = max(0.0, δf*(p[t]*(C[t] + ΔIN[t] + Igf[t]) + pk[t]*I[t] - Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] - i*Lf[t])) #   行動方程式
    #   企業利潤の配分
    Π[t] = p[t]*(C[t] + ΔIN[t] + Igf[t]) + pk[t]*I[t] - Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] - i*Lf[t] - Tf[t]   #   垂直一貫性の会計恒等式
    if Π[t] >= 0
        Πf[t] = ϵ*Π[t]  #   行動方程式
    else
        Πf[t] = Π[t]  #   行動方程式
    end
    Πi[t] = (Π[t] - Πf[t])*Ei[t]/(Ei[t] + Eb[t])    #   行動方程式
    Πb[t] = Π[t] - Πf[t] - Πi[t]    #   水平一貫性の会計恒等式
    #   金融機関の賃金の流動的な決定（モデルを簡単にするために、金融機関の収支を無理やり0にしようとした、不自然な形。）
    WbNb[t] = Πb[t] + ig*GBb[t] + i*L[t] #   垂直一貫性の会計恒等式
    WN[t] = Wf[t]*(Nc[t] + Nk[t]) + WbNb[t] + Wg[t]*Ng[t]   #   水平一貫性の会計恒等式
    #   可処分所得と確定した原価の計算
    YDw[t] = WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t]  #  定義式
    YDi[t] = Πi[t] + ig*GBi[t] - Tii[t] - Tei[t]   #  定義式
    if S[t] + Igf[t] + IN[t] == 0.0
        UC[t] = UC[tm1]
    else
        UC[t] = (UC[tm1]*IN[t] + (Wf[t]*(Nc[t] + Nk[t]) + Tv[t] + Tf[t] + Tef[t]))/(S[t] + Igf[t] + IN[t])   #   定義式
    end
    #   政府の貯蓄の計算
    GS[t] = -p[t]*Cg[t] + p[t]*Igg[t] - Wg[t]*Ng[t] - SS[t] + Tv[t] + Ti[t] + Te[t] + Tf[t] - ig*GB[t] #   垂直一貫性の会計恒等式
    #   純貸出と純資産の更新。
    NLw[t] = -p[t]*Cw[t] + WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t]    #   定義式
    NLi[t] = -p[t]*Ci[t] - Tii[t] - Tei[t] + Πi[t] + ig*GBi[t]  #   定義式
    NLf[t] = -p[t]*ΔIN[t] - pk[t]*I[t] + Πf[t]  #   定義式
    NLg[t] = -p[t]*Ig[t] + GS[t]   #   定義式
    tmp = abs(NLw[t+1] + NLi[t+1] + NLf[t+1] + NLg[t+1])
    if tmp > 0.00001
        println("会計的一貫性が崩れている p0, t=",t)
    end
    #   労働者の資金調達
    Lw[t+1] = η*WN[t]   #  行動方程式
    ΔLw[t] = Lw[t+1] - Lw[t]    #   ストックとフローの整合性の会計恒等式
    ΔMw[t] = NLw[t] + ΔLw[t]    #   垂直一貫性の会計恒等式
    Mw[t+1] = Mw[t] + ΔMw[t]    #   ストックとフローの整合性の会計恒等式
    #   企業の資金調達
    ΔMf[t] = NLf[t] #   仮置きの、垂直一貫性の会計恒等式。下の値の書き換えの過程で、ΔMf[t] - pe*ΔE[t] - ΔLf[t] = NLf[t] を常に満たすようにできればOK
    if Mf[t] > Lf[t]    #   返済
        ΔLf[t] -= 0.1*Lf[t] #   行動方程式
        ΔMf[t] -= 0.1*Lf[t] #   行動方程式
    end
    if Mf[t] + ΔMf[t] < Wf[t]*(Nc[t] + Nk[t])  #   流動資産が不足気味の時
        if I[t] > 0     #   資本が増加しているときは株式発行で対応
            ΔE[t] += pk[t]*I[t]/pe  #   行動方程式
            ΔMf[t] += pk[t]*I[t]  #   #   垂直一貫性の会計恒等式を維持するための操作(株買ったら手持ちの預金が同額減るともいう)
        end
        if ΔIN[t] > 0   #   在庫増加しているときは借入で対応
            ΔLf[t] += p[t]*ΔIN[t]  #   行動方程式
            ΔMf[t] += p[t]*ΔIN[t]    #   垂直一貫性の会計恒等式を維持するための操作(金借りたら同額預金が増えるとも言う)
        end
        if Mf[t] + ΔMf[t] < Wf[t]*(Nc[t] + Nk[t])  #   それでも流動資産が不足気味の時(シミュレーション開始時と、緊急用の分岐。この分岐はできるだけ使いたくない)
            short = Wf[t]*(Nc[t] + Nk[t]) - (Mf[t] + ΔMf[t])
            ΔE[t] += 0.5*short/pe  #   行動方程式
            ΔLf[t] += 0.5*short  #   行動方程式
            ΔMf[t] += short #   垂直一貫性の会計恒等式を維持するための操作
        end
    end
    Lf[t+1] = Lf[t] + ΔLf[t]    #   ストックとフローの整合性の会計恒等式
    E[t+1] = E[t] + ΔE[t]       #   ストックとフローの整合性の会計恒等式
    Mf[t+1] = Mf[t] + ΔMf[t]    #   ストックとフローの整合性の会計恒等式
    if abs(NLf[t] - ΔMf[t] + ΔLf[t] + pe*ΔE[t]) > 0.00001
        println("会計的一貫性が壊れている p1 t=",t)
    end
    #   投資家のポートフォリオ配分
    V, Y = NWi[t] + NLi[t], NLi[t]
    #   λ01 + λ11 + λ21 = 0, λ00 + λ10 + λ20 = 1　となるλを選んでいるため、垂直一貫性の会計恒等式が保証されている
    Mi[t+1] = λ00*V + λ01*Y     #   行動方程式
    Ei[t+1] = (λ10*V + λ11*Y)/pe    #   行動方程式
    GBi[t+1] = λ20*V + λ21*Y    #   行動方程式
    if GBi[t+1] < 0.0   #   GBi[t+1] を<0にさせないための緊急用の措置
        Mi[t+1] += -GBi[t+1]  #   行動方程式
        GBi[t+1] = 0.0  #   垂直統合性を維持するための式
    end
    if Ei[t+1] < 0.0
        Mi[t+1] += -pe*Ei[t+1]  #   行動方程式
        Ei[t+1] = 0.0   #   垂直統合性を維持するための式
    end
    if Ei[t+1] > E[t+1]
        Mi[t+1] += pe*(Ei[t+1] - E[t+1])
        Ei[t+1] = E[t+1]
    end
    if Mi[t+1] < 0.0
        Li[t+1] = -Mi[t+1]  #   行動方程式
        Mi[t+1] = 0.0       #   垂直統合性を維持するための式
    end
    ΔLi[t] = Li[t+1] - Li[t]    #   ストックとフローの整合性の会計恒等式
    ΔMi[t] = Mi[t+1] - Mi[t]    #   ストックとフローの整合性の会計恒等式
    ΔEi[t] = Ei[t+1] - Ei[t]    #   ストックとフローの整合性の会計恒等式
    ΔGBi[t] = GBi[t+1] - GBi[t] #   ストックとフローの整合性の会計恒等式
    if abs(NLi[t] - ΔMi[t] + ΔLi[t] - pe*ΔEi[t] - ΔGBi[t]) > 0.00001
        println("会計的一貫性が壊れている p2 t=",t)
    end
    #   金融機関のポートフォリオ。投資家の意思決定の結果を押し付けられるモデルにした。
    ΔEb[t] = ΔE[t] - ΔEi[t]     #   水平一貫性の会計恒等式
    Eb[t+1] = Eb[t] + ΔEb[t]  #   ストックとフローの整合性の会計恒等式
    ΔM[t] = ΔMw[t] + ΔMi[t] + ΔMf[t]    #   水平一貫性の会計恒等式
    M[t+1] = M[t] + ΔM[t]   #   ストックとフローの整合性の会計恒等式
    H[t+1] = ζ*M[t+1]   #   行動方程式　準備金制度
    ΔH[t] = H[t+1] - H[t]   #   ストックとフローの整合性の会計恒等式
    ΔL[t] = ΔLw[t] + ΔLi[t] + ΔLf[t] #   水平一貫性の会計恒等式
    L[t+1] = L[t] + ΔL[t]   #   ストックとフローの整合性の会計恒等式
    ΔGBb[t] = ΔM[t] - ΔL[t] - ΔH[t] - pe*ΔEb[t]   #   垂直一貫性の会計恒等式
    GBb[t+1] = GBb[t] + ΔGBb[t] #   ストックとフローの整合性の会計恒等式
    if abs(ΔM[t] - ΔL[t] - ΔH[t] - pe*ΔEb[t] - ΔGBb[t]) > 0.00001
        println("会計的一貫性が壊れている p3 t=",t)
    end
    #   統合政府のポートフォリオ配分。民間金融機関の意思決定の結果から内生的にHPMや国債を提供する
    ΔGB[t] = ΔGBi[t] + ΔGBb[t]    #   水平一貫性の会計恒等式
    GB[t+1] = GB[t] + ΔGB[t]    #   ストックとフローの整合性の会計恒等式
    #   隠された恒等式による、会計的一貫性の確認
    if abs(Cw[t] + Ci[t] + Cg[t] - C[t]) > 0.00001
        println("水平統合性が損なわれている p4 t=",t," Δ=",Cw[t] + Ci[t] + Cg[t] - C[t])
    end
    if abs(Igf[t] + Igg[t] - Ig[t]) > 0.00001
        println("水平統合性が損なわれている p5 t=",t," Δ=",Igf[t] + Igg[t] - Ig[t])
    end
    if abs(Wf[t]*(Nc[t]+Nk[t]) + WbNb[t] + Wg[t]*Ng[t] - WN[t]) > 0.00001
        println("水平統合性が損なわれている p6 t=",t," Δ=",Wf[t]*(Nc[t]+Nk[t]) + WbNb[t] + Wg[t]*Ng[t] - WN[t])
    end
    if abs(Tiw[t] + Tii[t] - Ti[t]) > 0.00001
        println("水平統合性が損なわれている p7 t=",t," Δ=",Tiw[t] + Tii[t] - Ti[t])
    end
    if abs(Tew[t] + Tei[t] + Tef[t] - Te[t]) > 0.00001
        println("水平統合性が損なわれている p8 t=",t," Δ=",Tew[t] + Tei[t] + Tef[t] - Te[t])
    end
    if abs(Πi[t] + Πf[t] + Πb[t] - Π[t]) > 0.00001
        println("水平統合性が損なわれている p9 t=",t," Δ=",Πi[t] + Πf[t] + Πb[t] - Π[t])
    end
    if abs(ΔMw[t] + ΔMi[t] + ΔMf[t] - ΔM[t]) > 0.00001
        println("水平統合性が損なわれている p10 t=",t," Δ=",ΔMw[t] + ΔMi[t] + ΔMf[t] - ΔM[t])
    end
    if abs(ΔLw[t] + ΔLi[t] + ΔLf[t] - ΔL[t]) > 0.00001
        println("水平統合性が損なわれている p11 t=",t," Δ=",ΔLw[t] + ΔLi[t] + ΔLf[t] - ΔL[t])
    end
    if abs(ΔEi[t] + ΔEb[t] - ΔE[t]) > 0.00001
        println("水平統合性が損なわれている p12 t=",t," Δ=",ΔEi[t] + ΔEb[t] - ΔE[t])
    end
    if abs(ΔGBi[t] + ΔGBb[t] - ΔGB[t]) > 0.00001
        println("水平統合性が損なわれている p13 t=",t," Δ=",ΔGBi[t] + ΔGBb[t] - ΔGB[t])
    end
    if abs(-p[t]*Cw[t] + WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t] - ΔMw[t] + ΔLw[t]) > 0.00001
        println("垂直統合性が損なわれている p14 t=",t," Δ=",-p[t]*Cw[t] + WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t] - ΔMw[t] + ΔLw[t])
    end
    if abs(-p[t]*Ci[t]- Tii[t] - Tei[t] + Πi[t] + ig*GBi[t] - ΔMi[t] + ΔLi[t] - pe*ΔEi[t] - ΔGBi[t]) > 0.00001
        println("垂直統合性が損なわれている p15 t=",t," Δ=",-p[t]*Ci[t]- Tii[t] - Tei[t] + Πi[t] + ig*GBi[t] - ΔMi[t] + ΔLi[t] - pe*ΔEi[t] - ΔGBi[t])
    end
    if abs(p[t]*C[t] + p[t]*ΔIN[t] + pk[t]*I[t] + p[t]*Igf[t] - Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] -Tf[t] - Π[t] - i*Lf[t]) > 0.00001
        println("垂直統合性が損なわれている p16 t=",t," Δ=",p[t]*C[t] + p[t]*ΔIN[t] + pk[t]*I[t] + p[t]*Igf[t] - Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] -Tf[t] - Π[t] - i*Lf[t])
    end
    if abs(-p[t]*ΔIN[t] - pk[t]*I[t] + Πf[t] - ΔMf[t] + ΔLf[t] + pe*ΔE[t]) > 0.00001
        println("垂直統合性が損なわれている p17 t=",t," Δ=",-p[t]*ΔIN[t] - pk[t]*I[t] + Πf[t] - ΔMf[t] + ΔLf[t] + pe*ΔE[t])
    end
    if abs(-WbNb[t] + Πb[t] + ig*GBb[t] + i*L[t]) > 0.00001
        println("垂直統合性が損なわれている p18 t=",t," Δ=",-WbNb[t] + Πb[t] + ig*GBb[t] + i*L[t])
    end
    if abs(ΔM[t] - ΔL[t] - ΔH[t] - pe*ΔEb[t] - ΔGBb[t]) > 0.00001
        println("垂直統合性が損なわれている p19 t=",t," Δ=",ΔM[t] - ΔL[t] - ΔH[t] - pe*ΔEb[t] - ΔGBb[t])
    end
    if abs(-p[t]*Cg[t] + p[t]*Igg[t] - Wg[t]*Ng[t] - SS[t] + Tv[t] + Ti[t] + Te[t] + Tf[t] - ig*GB[t] - GS[t]) > 0.00001
        println("垂直統合性が損なわれている p20 t=",t," Δ=",-p[t]*Cg[t] + p[t]*Igg[t] - Wg[t]*Ng[t] - SS[t] + Tv[t] + Ti[t] + Te[t] + Tf[t] - ig*GB[t] - GS[t])
    end
    if abs(NLg[t] + ΔH[t] + ΔGB[t]) > 0.00001#   隠された会計恒等式
        println("垂直統合性が損なわれている p21 t=",t," Δ=",NLg[t] + ΔH[t] + ΔGB[t])
    end
    #   BSMの会計恒等式を追加する必要があるかどうか未確認
    NWw[t+1] = Mw[t+1] - Lw[t+1]  #   ストックとフローの整合性の会計恒等式
    NWi[t+1] = Mi[t+1] + pe*Ei[t+1] + GBi[t+1] - Li[t+1]  #   ストックとフローの整合性の会計恒等式
    NWf[t+1] = pk[t+1]*K[t+1] + p[t+1]*IN[t+1] + Mf[t+1] - Lf[t+1] - pe*E[t+1]  #   ストックとフローの整合性の会計恒等式
    NWg[t+1] = p[t+1]*Kg[t+1] - H[t+1] - GB[t+1]  #   ストックとフローの整合性の会計恒等式
end

#   パラメータの変更・外生変数の変更


#   影響のシミュレーション


#   可視化
function plots_per_nation()
    plot(p[end-100:end-1].*Cg[end-100:end-1], label="Cg")
    plot!(p[end-100:end-1].*Igf[end-100:end-1], label="Igf")
    plot!(SS[end-100:end-1], label="SS")
    plot!(Wg[end-100:end-1].*Ng[end-100:end-1], label="Wg")
    savefig("figs/per_nation/GovExp.png")
    plot(Tv[end-100:end-1], label="Tv")
    plot!(Ti[end-100:end-1], label="Ti")
    plot!(Te[end-100:end-1], label="Te")
    plot!(Tf[end-100:end-1], label="Tf")
    savefig("figs/per_nation/Tax.png")
    plot(p[end-100:end-1].*C[end-100:end-1], label="C")
    plot!(p[end-100:end-1].*Cw[end-100:end-1], label="Cw")
    plot!(p[end-100:end-1].*Ci[end-100:end-1], label="Ci")
    plot!(p[end-100:end-1].*Cg[end-100:end-1], label="Cg")
    plot!(p[end-100:end-1].*S[end-100:end-1], label="S")
    plot!(p[end-100:end-1].*IN[end-100:end-1], label="IN")
    savefig("figs/per_nation/C_and_S.png")
    plot(NWw[end-100:end-1], label="NWw")
    plot!(NWi[end-100:end-1], label="NWi")
    plot!(NWf[end-100:end-1], label="NWf")
    plot!(NWg[end-100:end-1], label="NWg")
    savefig("figs/per_nation/NW.png")
    plot(NLw[end-100:end-1], label="NLw")
    plot!(NLi[end-100:end-1], label="NLi")
    plot!(NLf[end-100:end-1], label="NLf")
    plot!(NLg[end-100:end-1], label="NLg")
    savefig("figs/per_nation/NL.png")
    plot(Nk[end-100:end-1], label="Nk")
    plot!(Nc[end-100:end-1], label="Nc")
    plot!(Ng[end-100:end-1], label="Ng")
    savefig("figs/per_nation/N.png")
    plot(uk[end-100:end-1], label="uk")
    plot!(uc[end-100:end-1], label="uc")
    plot!(u[end-100:end-1], label="u")
    savefig("figs/per_nation/u.png")
    plot(ugk[end-100:end-1], label="ugk")
    plot!(ugc[end-100:end-1], label="ugc")
    plot!(ug[end-100:end-1], label="ug")
    savefig("figs/per_nation/ug.png")
    plot(I[end-100:end-1], label="I")
    plot!(Ig[end-100:end-1], label="Ig")
    plot!(Kg[end-100:end-1], label="Kg")
    plot!(K[end-100:end-1], label="K")
    savefig("figs/per_nation/K.png")
    plot(M[end-100:end-1], label="M")
    plot!(Mw[end-100:end-1], label="Mw")
    plot!(Mi[end-100:end-1], label="Mi")
    plot!(Mf[end-100:end-1], label="Mf")
    savefig("figs/per_nation/M.png")
    plot(Wf[end-100:end-1], label="Wf")
    plot!(Wg[end-100:end-1], label="Wg")
    savefig("figs/per_nation/W.png")
    plot(p[end-100:end-1], label="p")
    plot!(pk[end-100:end-1], label="pk")
    plot!(UC[end-100:end-1], label="UC")
    savefig("figs/per_nation/p.png")
    plot(YDw[end-100:end-1], label="YDw")
    plot!(YDi[end-100:end-1], label="YDi")
    savefig("figs/per_nation/YD.png")
    plot(Tv[end-100:end-1], label="Tv")
    plot!(Tf[end-100:end-1], label="Tf")
    savefig("figs/per_nation/Tv_and_Tf.png")
    plot(Tiw[end-100:end-1], label="Tiw")
    plot!(Tii[end-100:end-1], label="Tii")
    plot!(Ti[end-100:end-1], label="Ti")
    savefig("figs/per_nation/Ti.png")
    plot(Tew[end-100:end-1], label="Tew")
    plot!(Tei[end-100:end-1], label="Tei")
    plot!(Tef[end-100:end-1], label="Tef")
    plot!(Te[end-100:end-1], label="Te")
    savefig("figs/per_nation/Te.png")
    plot(C_demand[end-100:end-1], label="C_demand")
    plot!(C[end-100:end-1], label="C")
    plot!(Ce[end-100:end-1], label="Ce")
    savefig("figs/per_nation/C_demand")
    plot(Π[end-100:end-1], label="Π")
    plot!(Πi[end-100:end-1], label="Πi")
    plot!(Πf[end-100:end-1], label="Πf")
    plot!(Πb[end-100:end-1], label="Πb")
    savefig("figs/per_nation/Π")
    plot(E[end-100:end-1], label="E")
    plot!(Ei[end-100:end-1], label="Ei")
    plot!(Eb[end-100:end-1], label="Eb")
    savefig("figs/per_nation/E")
    plot(L[end-100:end-1], label="L")
    plot!(Lw[end-100:end-1], label="Lw")
    plot!(Li[end-100:end-1], label="Li")
    plot!(Lf[end-100:end-1], label="Lf")
    savefig("figs/per_nation/L")
end
function plots_per_worker()
    plot(SS[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="SS")
    savefig("figs/per_worker/SS.png")
    plot(Tv[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tv")
    plot!(Ti[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ti")
    plot!(Te[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Te")
    plot!(Tf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tf")
    savefig("figs/per_worker/Tax.png")
    plot(p[end-100:end-1].*C[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="C")
    plot!(p[end-100:end-1].*Cw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Cw")
    plot!(p[end-100:end-1].*Ci[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ci")
    plot!(p[end-100:end-1].*Cg[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Cg")
    plot!(p[end-100:end-1].*S[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="S")
    plot!(p[end-100:end-1].*IN[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="IN")
    savefig("figs/per_worker/C_and_S.png")
    plot(NWw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NWw")
    plot!(NWi[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NWi")
    plot!(NWf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NWf")
    plot!(NWg[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NWg")
    savefig("figs/per_worker/NW.png")
    plot(NLw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NLw")
    plot!(NLi[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NLi")
    plot!(NLf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NLf")
    plot!(NLg[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="NLg")
    savefig("figs/per_worker/NL.png")
    plot(I[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="I")
    plot!(Ig[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ig")
    plot!(Kg[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Kg")
    plot!(K[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="K")
    savefig("figs/per_worker/K.png")
    plot(M[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="M")
    plot!(Mw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Mw")
    plot!(Mi[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Mi")
    plot!(Mf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Mf")
    savefig("figs/per_worker/M.png")
    plot(Tv[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tv")
    plot!(Tf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tf")
    savefig("figs/per_worker/Tv_and_Tf.png")
    plot(Tiw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tiw")
    plot!(Tii[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tii")
    plot!(Ti[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ti")
    savefig("figs/per_worker/Ti.png")
    plot(Tew[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tew")
    plot!(Tei[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tei")
    plot!(Tef[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Tef")
    plot!(Te[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Te")
    savefig("figs/per_worker/Te.png")
    plot(C_demand[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="C_demand")
    plot!(C[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="C")
    plot!(Ce[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ce")
    savefig("figs/per_worker/C_demand")
    plot(Π[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Π")
    plot!(Πi[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Πi")
    plot!(Πf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Πf")
    plot!(Πb[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Πb")
    savefig("figs/per_worker/Π")
    plot(E[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="E")
    plot!(Ei[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Ei")
    plot!(Eb[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Eb")
    savefig("figs/per_worker/E")
    plot(L[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="L")
    plot!(Lw[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Lw")
    plot!(Li[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Li")
    plot!(Lf[end-100:end-1]./(Ng[end-100:end-1] + Nc[end-100:end-1] + Nk[end-100:end-1]), label="Lf")
    savefig("figs/per_worker/L")
end
plots_per_nation()
plots_per_worker()