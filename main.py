# coding=utf-8
import math
import pandas as pd


def cal_pmv(ta=None, tr=None, rh=None, v=None, clo=None, met=None, top=None, w=0):
    """
        This function is used to calculated Predicted Mean Vote(PMV) through environmental parameters
        cal_pmv(ta=None, tr=None, rh=None, v=None, clo=None, met=None, top=None, w=0)
        ta: Air temperature, °C
        tr: Mean radiation temperature, °C
        rh: Relative humidity, %
        v:  Relative air velocity, m/s
        clo: Clothing insulation, clo
        met: Metabolic rate, met
        top: operative temperature, °C
        When top is NOT None, it will replace values of ta AND tr
        w: External work, W/m2

        example 1:
        Out[1]: cal_pmv(ta=25, tr=25, rh=50, v=0.1, clo=0.5, met=1.0)
        Out[1]: -0.4025360559202205

        example 2:
        Out[1]: cal_pmv(to=25, rh=50, v=0.1, clo=0.5, met=1.0)
        Out[1]: -0.4025360559202205
    """
    if top is not None:
        ta = top
        tr = top
    m = 58.15 * met
    if clo < 0.5:
        fcl = 1 + 0.2 * clo
    else:
        fcl = 1.05 + .1 * clo
    icl = 0.155 * clo
    tsk = 35.7 - 0.028 * m
    pa = rh / 100 * math.exp(16.6536 - 4030.183 / (ta + 235))
    tcl = min(ta, tsk)
    hcf = 12.1 * v ** 0.5
    max_iteration = int(1e5)
    hc = hcf
    hr = 4
    for i in range(0, max_iteration):
        hcn = 2.38 * abs(ta - tcl) ** 0.25
        hc = max(hcf, hcn)
        tcl_k = tcl + 273
        tr_k = tr + 273
        hr = 3.96e-8 * (tcl_k ** 2 + tr_k ** 2) * (tcl_k + tr_k)
        tcl1 = (icl * fcl * (hc * ta + hr * tr) + tsk) / (1 + icl * fcl * (hr + hc))
        if abs(tcl - tcl1) < 1e-4:
            break
        tcl = tcl1

    r = fcl * hr * (tcl - tr)
    c = fcl * hc * (tcl - ta)
    c_res = 0.0014 * m * (34 - ta)
    e_res = 0.017 * m * (5.867 - pa)
    e_dif = 3.05 * (5.733 - 0.00699 * (m - w) - pa)
    e_rsw = max(0, 0.42 * (m - w - 58.15))
    tl = m - w - c - r - e_rsw - e_res - e_dif - c_res
    pmv = (0.303 * math.exp(-0.036 * m) + 0.028) * tl
    return pmv


def cal_ppd(pmv):
    """
            This function is used to calculated Predicted Percentage of Dissatisfaction(PPD) through PMV
            PPD = 100 - 95 * exp(-(0.03353 * PMV^4 + 0.2179 * PMV^2))
    """
    return 100 - 95 * math.exp(-(0.03353 * pmv ** 4 + 0.2179 * pmv ** 2))


class TwoNodeModel:
    met_factor = 58.2
    c_sw = 170
    c_dil = 120
    c_str = 0.5
    t_sk_0 = 33.7
    t_cr_0 = 36.8
    t_b_0 = 36.49
    skin_blood_flow_0 = 6.3
    sbc = 5.6697e-8
    alpha_0 = 0.1
    t_sk = t_sk_0
    t_cr = t_cr_0
    t_b = t_b_0
    skin_blood_flow = skin_blood_flow_0
    alpha = alpha_0
    h_sk = 0
    w = 0
    ps_sk = 0
    hr = 4.7
    met = 0
    work = 0

    def __init__(self, height=1.8, weight=70):
        self.height = height
        self.weight = weight
        self.a_du = 0.20247 * self.height ** 0.725 * self.weight ** 0.425

    def two_node_set(self, t_sk=t_sk_0, t_cr=t_cr_0,
                     skin_blood_flow=skin_blood_flow_0, alpha=alpha_0):
        self.t_sk = t_sk
        self.t_cr = t_cr
        self.skin_blood_flow = skin_blood_flow
        self.alpha = alpha

    def two_node_cal(self, ta=None, tr=None, v=None, rh=None, met=None, clo=None, b=101.325, top=None, icl=0.45,
                     work=0, duration_time=100, ashrae=False):
        if ashrae is True:
            self.weight = 69.9
            self.a_du = 1.8258
            duration_time = 60
        if top is not None:
            ta = top
            tr = top
        b = b*0.009869
        pw = rh/100*cal_pws(ta)
        v = max(v, 0.1)
        lr = 2.2/b
        rcl = clo*0.155
        fcl = 1 + 0.15*clo       # 也不知道为什么这里是0.15而后面要用0.25，未解之谜
        self.met = met
        self.work = work
        rm = met*self.met_factor
        m = met*self.met_factor
        if clo <= 0:
            w_crit = 0.38*v**(-0.29)
            icl = 1
        else:
            w_crit = 0.59*v**(-0.08)
        hc_n = 3*b**0.53
        hc_c = 8.600001*(v*b)**0.53
        hc = max(hc_n, hc_c)
        ht = hc + self.hr
        ra = 1/(fcl*ht)   # 空气层热阻
        top = (self.hr*tr + hc*ta)/ht
        tcl = top + (self.t_sk - top)/(ht*(ra+rcl))  # (tcl-top)*ht = (t_sk - top)/(ra+rcl) 似乎不太对
        e_sk = met * 0.1
        dry = 0
        p_wet = 0
        flag = 1
        tcl1 = tcl
        for i in range(0, duration_time):
            # 传热计算
            while 1:    # 迭代确定辐射换热系数
                if flag == 1:
                    tcl1 = tcl
                    self.hr = 4 * self.sbc * ((tcl + tr) / 2 + 273.15) ** 3 * 0.72
                    ht = hc + self.hr
                    ra = 1 / (fcl * ht)
                    top = (self.hr * tr + hc * ta) / ht
                tcl = (ra*self.t_sk + rcl*top)/(ra + rcl)    # (tcl-top)/ra = (t_sk-tcl)/rcl
                flag = 1
                if abs(tcl-tcl1) <= 1e-2:
                    break
            flag = 0
            dry = (self.t_sk - top)/(ra+rcl)
            hfcs = (self.t_cr - self.t_sk)*(5.28 + 1.163*self.skin_blood_flow)     # 血流换热量
            e_res = 0.0023*m*(44-pw)
            c_res = 0.0014*m*(34-ta)
            s_cr = m - hfcs - e_res - c_res - work
            s_sk = hfcs - dry - e_sk
            tc_sk = 0.97*self.alpha*self.weight
            tc_cr = 0.97*(1-self.alpha)*self.weight
            dt_sk = (s_sk*self.a_du)/(tc_sk*60)
            dt_cr = (s_cr*self.a_du)/(tc_cr*60)
            self.t_sk = self.t_sk + dt_sk
            self.t_cr = self.t_cr + dt_cr
            self.t_b = self.alpha*self.t_sk + (1 - self.alpha)*self.t_cr
            # 控制信号计算
            sk_sig = self.t_sk - self.t_sk_0
            warm_s = max(0.0, sk_sig)
            cold_s = max(0.0, -sk_sig)
            cr_sig = self.t_cr - self.t_cr_0
            warm_c = max(0.0, cr_sig)
            cold_c = max(0.0, -cr_sig)
            bd_sig = self.t_b-self.t_b_0
            warm_b = max(0.0, bd_sig)
            self.skin_blood_flow = (self.skin_blood_flow_0 + self.c_dil*warm_c)/(1 + self.c_str*cold_s)
            self.skin_blood_flow = max(0.5, min(90.0, self.skin_blood_flow))
            reg_sw = self.c_sw*warm_b*math.exp(warm_s/10.7)
            reg_sw = min(reg_sw, 500.0)
            e_rsw = 0.68*reg_sw
            rea = 1/(lr*fcl*hc)
            recl = rcl/(lr*icl)
            e_max = (cal_pws(self.t_sk) - pw)/(rea + recl)
            p_rsw = e_rsw/e_max
            p_wet = 0.06 + 0.94*p_rsw
            e_dif = p_wet*e_max - e_rsw
            if p_wet > w_crit:
                p_wet = w_crit
                p_rsw = w_crit/0.94
                e_rsw = p_rsw*e_max
                e_dif = 0.06*(1-p_rsw)*e_max
            if e_max < 0:
                e_dif = 0
                e_rsw = 0
                p_wet = w_crit
            e_sk = e_rsw + e_dif
            m_shiv = 19.4*cold_s*cold_c
            m = rm + m_shiv
            self.alpha = 0.0417737+0.7451833/(self.skin_blood_flow+0.585417)
        self.h_sk = dry + e_sk
        self.w = p_wet
        self.ps_sk = cal_pws(self.t_sk)

    def get_set(self,  b=101.325, work=0):
        b = b * 0.009869
        lr = 2.2 / b
        hr_s = self.hr
        if self.met < 0.85:
            hc_s = 3
        else:
            hc_s = max(5.66 * (self.met - 0.85) ** 0.39, 3)
        clo_s = 1.52 / ((self.met - work / self.met_factor) + 0.6944) - 0.1835
        ht_s = hc_s + hr_s
        rcl_s = 0.155 * clo_s
        facl_s = 1 + 0.25 * clo_s       #
        im_s = 0.45
        ra_s = 1 / (facl_s * ht_s)
        rea_s = 1 / (lr * facl_s * hc_s)
        fcl_s = 1 / (1 + 0.155 * facl_s * ht_s * clo_s)
        icl_s = im_s * hc_s / ht_s * (1 - fcl_s) / (hc_s / ht_s - fcl_s * im_s)
        recl_s = rcl_s / (lr * icl_s)
        hd_s = 1 / (ra_s + rcl_s)
        he_s = 1 / (rea_s + recl_s)
        # 用牛顿法求SET
        dt = 1e-4
        dx = 100
        set0 = self.t_sk - self.h_sk / hd_s
        while abs(dx) > 1e-3:
            err1 = self.h_sk - hd_s * (self.t_sk - set0) - self.w * he_s * (self.ps_sk - 0.5 * cal_pws(set0))
            err2 = self.h_sk - hd_s * (self.t_sk - (set0 + dt)) - self.w * he_s * (self.ps_sk - 0.5 * cal_pws(set0 + dt))
            dx = err1 / (err2 - err1) * dt
            set0 = set0 - dx
        return set0


def cal_pws(t):
    return math.exp(18.6686-4030.183/(t+235))


print("如果要计算PMV、PPD和SET\n数据组织形式参考template.csv文件，并且和本exe放在同一目录下")
while True:
    str = input("要继续请输入Y，退出请输入Q：")
    if str == "Q":
        break
    else:
        filename = input("请输入需要计算数据表格的文件名\t退出请输入Q\n（不包含扩展名，最好不要有中文，否则也不知道会不会有bug）：")
        path = "./" + filename + ".csv"
        output_path = "./" + filename + "_output.csv"

        f = pd.read_csv(path)
        tnm = TwoNodeModel()
        f["PMV"] = 0
        f["PPD"] = 0
        f["SET"] = 0

        for index, row in f.iterrows():
            ta = f.iloc[index, 0]
            tr = f.iloc[index, 1]
            rh = f.iloc[index, 2]
            v = f.iloc[index, 3]
            clo = f.iloc[index, 4]
            met = f.iloc[index, 5]
            tnm.two_node_set()
            tnm.two_node_cal(ta=ta, tr=tr, rh=rh, v=v, clo=clo, met=met, ashrae=True)
            pmv = cal_pmv(ta=ta, tr=tr, rh=rh, v=v, clo=clo, met=met)
            f.iloc[index, 6] = pmv
            f.iloc[index, 7] = cal_ppd(cal_pmv(ta=ta, tr=tr, rh=rh, v=v, clo=clo, met=met))
            f.iloc[index, 8] = tnm.get_set()

        f.to_csv(output_path, index=False)
        print("计算完成，好耶")


