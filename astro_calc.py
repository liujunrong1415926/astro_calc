import numpy as np
import astropy.constants as consts
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import chi2, norm
import streamlit as st

# =========================
# 常数
# =========================
au = consts.au.cgs.value
c = consts.c.cgs.value
day = 86400
G = consts.G.cgs.value
h_Gauss = consts.h.cgs.value
Jy = 1e-23
k_B = consts.k_B.cgs.value
L_sun = consts.L_sun.cgs.value
MeV = 1e6 * 1.6e-19 * 1e7
m_e = consts.m_e.cgs.value
micron = 1e-4
m_p = consts.m_p.cgs.value
Mpc = 1e6 * consts.pc.cgs.value
M_sun = consts.M_sun.cgs.value
M_earth = consts.M_earth.cgs.value
pc = consts.pc.cgs.value
R_sun = consts.R_sun.cgs.value
sigma_sb = consts.sigma_sb.cgs.value
sigma_T = consts.sigma_T.cgs.value
year = 365 * 86400

qe_SI = 1.602e-19
h_SI = 6.624e-34


# =========================
# 计算函数
# =========================
def z_to_dL(z):
    cosm1 = FlatLambdaCDM(H0=70, Om0=0.3)
    D_mpc = cosm1.luminosity_distance(z).value
    Dist = Mpc * D_mpc
    return {
        "z": z,
        "Dist_Mpc": D_mpc,
        "Dist_cm": Dist
    }


def flux_to_lum(flux=1e-13, Dist=None, z=None):
    if z is not None:
        Dist = z_to_dL(z)["Dist_cm"]
    if Dist is None:
        raise ValueError("必须提供 Dist 或 z 中的一个")
    L = flux * 4 * np.pi * Dist**2
    return {
        "flux": flux,
        "Dist_cm": Dist,
        "L": L
    }


def TS_calc(TS=None, sig=None, confidence=None, df=2):
    o = {"dof": df}

    if sig is not None:
        o["sig"] = float(sig)
        o["confidence"] = float(1 - 2 * (1 - norm.cdf(sig)))
        o["TS"] = float(chi2.ppf(o["confidence"], df=df))

    if TS is not None:
        o["TS"] = float(TS)
        o["confidence"] = float(chi2(df).cdf(TS))
        o["sig"] = float(norm.ppf(1 - (1 - o["confidence"]) / 2))

    if confidence is not None:
        o["confidence"] = float(confidence)
        o["TS"] = float(chi2.ppf(confidence, df=df))
        o["sig"] = float(norm.ppf(1 - (1 - confidence) / 2))

    return o


# =========================
# 页面设置
# =========================
st.set_page_config(page_title="Astro Calculator", layout="centered")
st.title("Astro Calculator")
# st.write("把你的 Python 计算封装成一个本地网页工具。")

tab1, tab2, tab3 = st.tabs([
    "z → dL",
    "flux → luminosity",
    "TS / sig / confidence"
])

# =========================
# 1. z_to_dL
# =========================
with tab1:
    st.subheader("红移转光度距离")
    z = st.number_input("输入红移 z", value=0.00254333, format="%.8f")
    if st.button("计算 dL"):
        out = z_to_dL(z)
        st.success("计算完成")
        st.write(f"z = {out['z']:.8f}")
        st.write(f"Dist = {out['Dist_Mpc']:.6e} Mpc")
        st.write(f"Dist = {out['Dist_cm']:.6e} cm")


# =========================
# 2. flux_to_lum
# =========================
with tab2:
    st.subheader("通量转光度")
    flux = st.number_input("输入 flux (erg s^-1 cm^-2)", value=1e-13, format="%.6e")

    mode = st.radio("选择输入方式", ["用距离 Dist", "用红移 z"])

    if mode == "用距离 Dist":
        Dist = st.number_input("输入 Dist (cm)", value=10.9 * Mpc, format="%.6e")
        if st.button("计算光度"):
            out = flux_to_lum(flux=flux, Dist=Dist)
            st.success("计算完成")
            st.write(f"flux = {out['flux']:.6e} erg s^-1 cm^-2")
            st.write(f"Dist = {out['Dist_cm']:.6e} cm")
            st.write(f"L = {out['L']:.6e} erg s^-1")

    else:
        z2 = st.number_input("输入红移 z", value=0.0025, format="%.8f")
        if st.button("按 z 计算光度"):
            out = flux_to_lum(flux=flux, z=z2)
            st.success("计算完成")
            st.write(f"flux = {out['flux']:.6e} erg s^-1 cm^-2")
            st.write(f"Dist = {out['Dist_cm']:.6e} cm")
            st.write(f"L = {out['L']:.6e} erg s^-1")


# =========================
# 3. TS_calc
# =========================
with tab3:
    st.subheader("TS / 显著性 / confidence 换算")

    input_type = st.selectbox(
        "选择已知量",
        ["TS", "sig", "confidence"]
    )

    df = st.number_input("自由度 df", min_value=1, value=2, step=1)

    if input_type == "TS":
        TS = st.number_input("输入 TS", value=54.0, format="%.6f")
        if st.button("由 TS 计算"):
            out = TS_calc(TS=TS, df=df)
            st.success("计算完成")
            st.write(f"sig = {out['sig']:.6f}")
            st.write(f"confidence = {out['confidence']:.10f}")
            st.write(f"1 - confidence = {1 - out['confidence']:.10e}")
            st.write(f"TS = {out['TS']:.6f}")
            st.write(f"dof = {out['dof']}")

    elif input_type == "sig":
        sig = st.number_input("输入 sig", value=4.0, format="%.6f")
        if st.button("由 sig 计算"):
            out = TS_calc(sig=sig, df=df)
            st.success("计算完成")
            st.write(f"sig = {out['sig']:.6f}")
            st.write(f"confidence = {out['confidence']:.10f}")
            st.write(f"1 - confidence = {1 - out['confidence']:.10e}")
            st.write(f"TS = {out['TS']:.6f}")
            st.write(f"dof = {out['dof']}")

    else:
        confidence = st.number_input("输入 confidence", value=0.95, format="%.10f")
        if st.button("由 confidence 计算"):
            out = TS_calc(confidence=confidence, df=df)
            st.success("计算完成")
            st.write(f"sig = {out['sig']:.6f}")
            st.write(f"confidence = {out['confidence']:.10f}")
            st.write(f"1 - confidence = {1 - out['confidence']:.10e}")
            st.write(f"TS = {out['TS']:.6f}")
            st.write(f"dof = {out['dof']}")


# =========================
# 页脚测试说明
# =========================
st.markdown("---")
st.markdown("你原来测试的例子：`TS=54, df=2` 和 `TS=71, df=2`，可直接在第三个标签页输入。")
