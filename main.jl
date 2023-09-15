#   関数定義
function btw(X::Float64, Y::Float64, Z::Float64)
    if X > Z    #   不正入力に対してはエラーを吐く値を返す
        return nothing
    end
    return max(X, min(Y, Z))
end

#   パラメータ
T = 100
ut = 0.8    #   目標稼働率
γ1, γ2, γ3 = 0.1, 0.1, 0.2

λ00, λ10, λ20, λ01, λ11, λ12 = 0.4, 0.4, 0.2, -1.0, 1.0, 0.0
λ21 = -(λ01 + λ11)
λ20 = 1 - (λ00 + λ10)

#   変数定義
m, p, pk = zeros(T), zeros(T), zeros(T)
Wf, Wg = zeros(T), zeros(T)
u, uc, uk = zeros(T), zeros(T), zeros(T)
ue = zeros(T)
UCe, UC = zeros(T), zeros(T)


#   初期値設定
m[1], p[1], pk[1] = 0.3, 1.0, 1.0
Wf[1], Wg[1] = 1.0, 1.0
u[1], uc[1], uk[1] = 0.8, 0.8, 0.8
ue[1] = 0.8
UCe[1], UC[1] = 1.0, 1.0
    #   ストックの変数の初期値は必ず会計的一貫性を持つように設定すること

#   定常状態に至るまでシミュレーション
for t = 2:T//2
    #   期待値計算
    ue[t] = λe*u[t-1] + (1 - λe)*ue[t-1]
    UCe[t] = λe*UC[t-1] + (1 - λe)*UCe[t-1]
    Ie[t] = λe*I[t-1] + (1 - λe)*Ie[t-1]
    Ke[t] = λe*K[t-1] + (1 - λe)*Ke[t-1]
    Kge[t] = λe*Kg[t-1] + (1 - λe)*Kge[t-1]
    Ce[t] = λe*C[t-1] + (1 - λe)*Ce[t-1]
    INe[t] = λe*IN[t-1] + (1 - λe)*INe[t-1]
    YDwe[t] = λe*YDw[t-1] + (1 - λe)*YDwe[t-1]
    YDie[t] = λe*YDi[t-1] + (1 - λe)*YDie[t-1]
    NWwe[t] = λe*NWw[t-1] + (1 - λe)*NWwe[t-1]
    NWie[t] = λe*NWi[t-1] + (1 - λe)*NWie[t-1]
    #   価格設定
    m[t] = m[t-1]*((1-γ2) + γ2*ue[t]/ut)    #   行動方程式
    UCe[t] = (UC[t-1]*IN[t]+(Wf[t]+Tve[t]+Tfe[t]+Tefe[t]))/(S[t]+IN[t]) #   行動方程式
    p[t] = (1+m[t])*UCe[t]  #   行動方程式
    pk[t] = Wf[t]*Nk[t]/Ie[t]   #   行動方程式
    Wf[t] = Wf[t-1]*((1-γ1) + γ1*ue[t]/ut)  #   行動方程式
    Wg[t] = Wg[t-1]*((1-γ3) + γ3*Wf[t-1]/Wg[t-1])   #   行動方程式
    #   雇用者数設定
    Nc[t] = btw((1-β4)*Nc[t-1], γc*β1*Ke[t], (1+β3)*Nc[t-1])    #   行動方程式
    Nk[t] = btw((1-β4)*Nk[t-1], γk*(1-β1)*Ke[t], (1+β3)*Nk[t-1])    #   行動方程式
    Ng[t] = btw((1-β4)*Ng[t-1], γg*βg3*Kge[t], (1+β3)*Ng[t-1])  #   行動方程式
    #   生産活動の意思決定
    Igf[t] = (β0+βg)*Kg[t-1]/(1+βg2)    #   行動方程式
    Igg[t] = (βg2*βg - β0)*Kg[t-1]/(1+βg2)  #   行動方程式
    Ig[t] = Igf[t] + Igg[t] #   水平一貫性の会計恒等式
    ug[t] = (Igg[t] + β0*Kg[t-1])/(βgk*βg3*Kg[t-1]) #   定義式
    S[t] = btw(0, Ce[t] - INe[t] + β2*Ce[t], btw(0, βc*Nc[t]/γc - Igf[t], βc*β1*K[t-1] - Igf[t]))   #   行動方程式
    I[t] = β0*K[t-1] + btw(0, ((S+Igf)/(βc*β1*βk) - ut)*K[t-1] + β0*K[t-1], βk*min((1-β1)*K[t-1],Nk[t]/γk, β5*Mf[t] + β6*(Π[t-1] - pk[t]*I[t-1])))  #   行動方程式
    uc[t] = (S[t] + Igf[t])/(βc*β1*K[t-1])  #   定義式
    uk[t] = (I[t] + β0*K[t-1])/(βk*(1-β1)*K[t-1])   #   定義式
    u[t] = (S[t] + Igf[t])/(βc*K[t-1]) + (I[t] + β0*K[t-1])/(βk*K[t-1]) #   定義式
    #   支出の意思決定
    Cw[t] = (α1*YDwe[t] + α2*NWwe[t])/p[t]  #   行動方程式
    Ci[t] = (α3*YDie[t] + α4*NWie[t])/p[t]  #   行動方程式
    Cg[t] = Cg0*(1 + α5)^t  #   行動方程式
    C[t] = Cw[t] + Ci[t] + Cg[t]    #   水平一貫性の会計恒等式
    SS[t] = SS0*(1 + α5)^t  #   行動方程式
    IN[t+1] = IN[t] + S[t] - C[t]   #   行動方程式
    #   税額の決定  （循環参照を避けるため、仕方なく、１期前の値から計算する。税額は後から確定するんだから問題ないってな言い訳もできるっちゃできる）
    Tv[t] = δv*(p[t-1]*C[t-1] + Wf[t-1]*(Nc[t-1] + Nk[t-1]) + pk[t-1]*Igf[t-1])   #   行動方程式
    Tiw[t] = δiw*WN[t-1]  #   行動方程式
    Tii[t] = δii*(Πi[t-1] + Di[t-1] + ig*GBi[t-1])    #   行動方程式
    Ti[t] = Tiw[t-1] + Tii[t-1] #   水平一貫性の会計恒等式
    Tew[t] = δew*max(0.0, Mw[t-1] - Lw[t-1])  #   行動方程式
    Tei[t] = δew*Mi + δei*(pe*Ei[t-1]+GBi[t-1]) #   行動方程式
    Tef[t] = δef*K[t-1] #   行動方程式
    Te[t] = Tew[t] + Tei[t] + Tef[t]    #   水平一貫性の会計恒等式
    Tf[t] = δf*(p[t]*C[t] + p[t]*ΔIN[t] + pk[t]*I[t] + pk[t]*Igf[t] + Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] - i*Lf[t]) #   行動方程式
    #   企業利潤の配分
    Π[t] = (1 - δf)*(p[t]*C[t] + p[t]*ΔIN[t] + pk[t]*I[t] + pk[t]*Igf[t] + Wf[t]*(Nc[t] + Nk[t]) - Tv[t] - Tef[t] - i*Lf[t])   #   垂直一貫性の会計恒等式
    Πf[t] = ϵ*Π[t]  #   行動方程式
    Πi[t] = (Π[t] - Πf[t])*Ei[t]/(Ei[t] + Eb[t])    #   行動方程式
    Πb[t] = Π[t] - Πf[t] - Πi[t]    #   水平一貫性の会計恒等式
    #   金融機関の賃金の流動的な決定（モデルを簡単にするために、金融機関の収支を無理やり0にしようとした、不自然な形。）
    WbNb[t] = Πb[t] + ig*GBb[t] + i*L[t] #   垂直一貫性の会計恒等式
    WN[t] = Wf[t]*(Nc[t] + Nk[t]) + WbNb[t] + Wg[t]*Ng[t]   #   水平一貫性の会計恒等式
    #   可処分所得の計算
    YDw[t] = WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t]  #  定義式
    YDi[t] = Πi + ig*GBi[t] - Tii[t] - Tei[t]   #  定義式
    #   政府の貯蓄の計算
    GS[t] = -p[t]*Cg[t] + pk[t]*Igg[t] - Wg[t]*Ng[t] - SS[t] + Tv[t] + Ti[t] + Te[t] + Tf[t] - ig*GB[t]
    #   純貸出と純資産の更新。
    NLw[t] = -p[t]*Cw[t] + WN[t] + SS[t] - Tiw[t] - Tew[t] - i*Lw[t]    #   定義式
    NLi[t] = -p[t]*Ci[t] - Tii[t] - Tei[t] + Πi[t] + ig*GBi[t]  #   定義式
    NLf[t] = -p[t]*ΔIN[t] - pk[t]*I[t] + Πf[t]  #   定義式
    NLg[t] = -pk[t]*Ig[t] + GS[t]   #   定義式
    NWw[t+1] = NWw[t] + NLw[t]  #   ストックとフローの整合性の会計恒等式(キャピタルゲインを考えなくてよいので成り立つ)
    NWi[t+1] = NWi[t] + NLi[t]  #   ストックとフローの整合性の会計恒等式(キャピタルゲインを考えなくてよいので成り立つ)
    NWf[t+1] = NWf[t] + NLf[t]  #   ストックとフローの整合性の会計恒等式(キャピタルゲインを考えなくてよいので成り立つ)
    NWg[t+1] = NWg[t] + NLg[t]  #   ストックとフローの整合性の会計恒等式(キャピタルゲインを考えなくてよいので成り立つ)
    #   労働者の資金調達
    Lw[t+1] = η*WN[t]   #  行動方程式
    ΔLw[t] = Lw[t+1] - Lw[t]    #   ストックとフローの整合性の会計恒等式
    ΔMw[t] = NLw[t] + ΔLw[t]    #   垂直一貫性の会計恒等式
    Mw[t+1] = Mw[t] + ΔMw[t]    #   ストックとフローの整合性の会計恒等式
    #   企業の資金調達
    if Mf[t] > Lf[t]    #   返済
        ΔLf[t] -= 0.1*Lf[t] #   行動方程式
    end
    ΔMf[t] = NLf[t] #   仮置きの、垂直一貫性の会計恒等式。下で、ΔMf[t] - pe*ΔE[t] - ΔLf[t] = NLf[t] を常に満たすようにできればOK
    if Mf[t+1] < Wf[t]*(Nc[t] + Nk[t])  #   流動資産が不足気味の時
        if ΔIN[t] > 0   #   在庫増加しているときは借入で対応
            ΔLf[t] += p*ΔIN[t]  #   行動方程式
            ΔMf[t] += ΔLf[t]    #   垂直一貫性の会計恒等式を維持するための操作(金借りたら同額預金が増えるとも言う)
        end
        if I[t] > 0     #   資本が増加しているときは株式発行で対応
            ΔE[t] = pk[t]*I[t]/pe  #   行動方程式
            ΔMf[t] += pe*ΔE[t]  #   #   垂直一貫性の会計恒等式を維持するための操作(株買ったら手持ちの預金が同額減るともいう)
        end
    end
    Lf[t+1] = Lf[t] + ΔLf[t]    #   ストックとフローの整合性の会計恒等式
    E[t+1] = E[t] + ΔE[t]       #   ストックとフローの整合性の会計恒等式
    Mf[t+1] = Mf[t] + ΔMf[t]    #   ストックとフローの整合性の会計恒等式
    #   投資家のポートフォリオ配分
    if NWi[t+1] <= 0.0  #   エラーを吐かせないための緊急用の措置。通常はこの条件で落ち着かないようにしなければならない
        Li[t+1] = NWi[t+1]  #   行動方程式
    else
        #   λ01 + λ11 + λ21 = 0, λ00 + λ10 + λ20 = 1　となるλを選んでいるため、垂直一貫性の会計恒等式が保証されている
        Mi[t+1] = λ00*NWi[t+1] + λ01*NLi[t]     #   行動方程式
        Ei[t+1] = (λ10*NWi[t+1] + λ11*NLi[t])/pe    #   行動方程式
        GBi[t+1] = λ20*NWi[t+1] + λ21*NLi[t]    #   行動方程式
        if min(Mi[t+1], Ei[t+1], GBi[t+1]) < 0   #   Mi[t+1], Ei[t+1], GBi[t+1] を<0にさせないための緊急用の措置
            Mi[t+1] = λ00*NWi[t+1]  #   行動方程式
            Ei[t+1] = λ10*NWi[t+1]  #   行動方程式
            GBi[t+1] = λ20*NWi[t+1] #   行動方程式
        end
    end
    ΔLi[t] = Li[t+1] - Li[t]    #   ストックとフローの整合性の会計恒等式
    ΔMi[t] = Mi[t+1] - Mi[t]    #   ストックとフローの整合性の会計恒等式
    ΔEi[t] = Ei[t+1] - Ei[t]    #   ストックとフローの整合性の会計恒等式
    ΔGBi[t] = GBi[t+1] - GBi[t] #   ストックとフローの整合性の会計恒等式
    #   金融機関のポートフォリオ。投資家の意思決定の結果を押し付けられるモデルにした。
    ΔEb[t] = ΔE[t] - ΔEi[t]     #   水平一貫性の会計恒等式
    Eb[t+1] = Eb[t] - ΔEb[t]  #   ストックとフローの整合性の会計恒等式
    ΔM[t] = ΔMw[t] + ΔMi[t] + ΔMf[t]    #   水平一貫性の会計恒等式
    M[t+1] = M[t] + ΔM[t]   #   ストックとフローの整合性の会計恒等式
    H[t+1] = ζ*M[t+1]   #   行動方程式　準備金制度
    ΔH[t] = H[t+1] - H[t]   #   ストックとフローの整合性の会計恒等式
    ΔGBb[t] = -ΔM[t] + ΔL[t] + ΔH[t] + pe*ΔEb[t]   #   垂直一貫性の会計恒等式
    GBb[t+1] = GBb[t] + ΔGBb[t] #   ストックとフローの整合性の会計恒等式
    ΔL[t] = ΔLw[t] + ΔLi[t] + ΔLf[t] #   水平一貫性の会計恒等式
    L[t+1] = L[t] + ΔL[t]   #   ストックとフローの整合性の会計恒等式
    #   統合政府のポートフォリオ配分。民間金融機関の意思決定の結果から内生的にHPMや国債を提供する
    ΔGB[t+1] = ΔGBi[t+1] + ΔGBb[t+1]    #   水平一貫性の会計恒等式
    #   隠された恒等式による、会計的一貫性の確認
    tmp = NLg[t] + ΔH[t] - ΔGB[t]   #   隠された会計恒等式
    if abs(tmp) > 0.001
        println("会計的一貫性が壊れている。t = $t")
    end
    
    #   BSMの会計恒等式を追加する必要があるかどうか未確認
end

#   パラメータの変更・外生変数の変更


#   影響のシミュレーション

#   可視化
    #   画像ファイルのアウトプットの方法を調べる