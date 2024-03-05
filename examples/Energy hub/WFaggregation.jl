using DelimitedFiles, SymEngine


s = symbols("s")

WT_Sbase = 5 # MW
WT_vACbase = 0.6
WF_vACbase = 66
WT_Zbase = WT_vACbase^2/WT_Sbase
Powf = 5
Qowf = 0

Powf_total = 50


net = @network begin

    WT_1 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_1 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_1 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_receiving = ac_source(V = WF_vACbase/sqrt(3), P = -Powf_total, P_min = -100, P_max = 100, Q_max = 100, Q_min = -100, pins = 3, transformation = true)

    G_sending_1[2.1] == gndd
    G_sending_1[2.2] == gndq

    G_receiving[2.1] == gndd
    G_receiving[2.2] == gndq

    G_sending_1[1.1] == MMC_sending_1[2.1]
    G_sending_1[1.2] == MMC_sending_1[2.2]

    MMC_sending_1[1.1] == WT_1[1.1]

    WT_2 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_2 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_2 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_2[2.1] == gndd
    G_sending_2[2.2] == gndq

    G_sending_2[1.1] == MMC_sending_2[2.1]
    G_sending_2[1.2] == MMC_sending_2[2.2]

    MMC_sending_2[1.1] == WT_2[1.1]

    WT_3 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_3 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_3 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_3[2.1] == gndd
    G_sending_3[2.2] == gndq

    G_sending_3[1.1] == MMC_sending_3[2.1]
    G_sending_3[1.2] == MMC_sending_3[2.2]

    MMC_sending_3[1.1] == WT_3[1.1]

    WT_4 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_4 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_4 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_4[2.1] == gndd
    G_sending_4[2.2] == gndq

    G_sending_4[1.1] == MMC_sending_4[2.1]
    G_sending_4[1.2] == MMC_sending_4[2.2]

    MMC_sending_4[1.1] == WT_4[1.1]

    WT_5 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_5 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_5 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_5[2.1] == gndd
    G_sending_5[2.2] == gndq

    G_sending_5[1.1] == MMC_sending_5[2.1]
    G_sending_5[1.2] == MMC_sending_5[2.2]

    MMC_sending_5[1.1] == WT_5[1.1]

    WT_6 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_6 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_6 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_6[2.1] == gndd
    G_sending_6[2.2] == gndq

    G_sending_6[1.1] == MMC_sending_6[2.1]
    G_sending_6[1.2] == MMC_sending_6[2.2]

    MMC_sending_6[1.1] == WT_6[1.1]

    WT_7 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_7 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_7 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_7[2.1] == gndd
    G_sending_7[2.2] == gndq

    G_sending_7[1.1] == MMC_sending_7[2.1]
    G_sending_7[1.2] == MMC_sending_7[2.2]

    MMC_sending_7[1.1] == WT_7[1.1]

    WT_8 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_8 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_8 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_8[2.1] == gndd
    G_sending_8[2.2] == gndq

    G_sending_8[1.1] == MMC_sending_8[2.1]
    G_sending_8[1.2] == MMC_sending_8[2.2]

    MMC_sending_8[1.1] == WT_8[1.1]

    WT_9 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_9 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_9 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_9[2.1] == gndd
    G_sending_9[2.2] == gndq

    G_sending_9[1.1] == MMC_sending_9[2.1]
    G_sending_9[1.2] == MMC_sending_9[2.2]

    MMC_sending_9[1.1] == WT_9[1.1]

    WT_10 = tlc(Vᵈᶜ = 1.3, Vₘ = WT_vACbase/sqrt(3), Lᵣ = 0.15*WT_Zbase/2/pi/50, Rᵣ = 0.02*WT_Zbase,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        P_max = 100, P_min = -100, P = Powf, Q = Qowf, Q_max = 100, Q_min = -100,
        occ = PI_control(Kₚ = 0.254647908947033, Kᵢ = 0.8),
        pll = PI_control(Kₚ = 0.795774715459477, Kᵢ = 31.830988618379067, ω_f = 2*pi*80),
        v_meas_filt = PI_control(ω_f = 1e4),
        i_meas_filt = PI_control(ω_f = 1e4),
        p = PI_control(Kₚ = 0.01, Kᵢ = 10),
        q = PI_control(Kₚ = 0.01, Kᵢ = 10)
        )

    # Extra AC/DC converter to represent DC voltage and get power flow converged 
    MMC_sending_10 = mmc(Vᵈᶜ = 1.3, vDCbase = 1.3, Vₘ = WT_vACbase/sqrt(3), 
        P_max = 100, P_min = -100, P = -Powf, Q = 0, Q_max = 100, Q_min = -100,
        vACbase_LL_RMS = WT_vACbase, Sbase = WT_Sbase,
        Lᵣ = 0.15 * WT_Zbase/2/pi/50, Rᵣ = 0.02 * WT_Zbase, Lₐᵣₘ = 0.0383 * WT_Zbase/2/pi/50, Rₐᵣₘ =0.026 *WT_Zbase, N = 10 , Cₐᵣₘ = 10e-3,
        occ = PI_control(Kₚ = 0.7691, Kᵢ = 522.7654),
        ccc = PI_control(Kₚ = 0.1048, Kᵢ = 48.1914),
        pll = PI_control(Kₚ = 0.28, Kᵢ = 12.5664),
        dc = PI_control(Kₚ = 5, Kᵢ = 15),
        q = PI_control(Kₚ = 0.1, Kᵢ = 31.4159)
        )

    G_sending_10 = ac_source(V = WT_vACbase/sqrt(3), P = Powf, P_min = -10, P_max = 10, Q_max = 10, Q_min = -10, pins = 3, transformation = true)

    G_sending_10[2.1] == gndd
    G_sending_10[2.2] == gndq

    G_sending_10[1.1] == MMC_sending_10[2.1]
    G_sending_10[1.2] == MMC_sending_10[2.2]

    MMC_sending_10[1.1] == WT_10[1.1]

    # Individual WT transformers
    t1 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t2 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t3 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t4 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t5 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t6 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t7 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t8 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t9 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)
    t10 = transformer(n = 0.69/33 , Lₚ = (0.15*(WT_vACbase^2/3)/(2*pi*50))/2, Rₚ = 0, Rₛ = 0, Lₛ = (0.15*(WF_vACbase^2/3)/(2*pi*50))/2, pins = 3, transformation = true)

    WT_1[2.1] == t1[1.1]
    WT_1[2.2] == t1[1.2]

    WT_2[2.1] == t2[1.1]
    WT_2[2.2] == t2[1.2]

    WT_3[2.1] == t3[1.1]
    WT_3[2.2] == t3[1.2]

    WT_4[2.1] == t4[1.1]
    WT_4[2.2] == t4[1.2]

    WT_5[2.1] == t5[1.1]
    WT_5[2.2] == t5[1.2]

    WT_6[2.1] == t6[1.1]
    WT_6[2.2] == t6[1.2]

    WT_7[2.1] == t7[1.1]
    WT_7[2.2] == t7[1.2]

    WT_8[2.1] == t8[1.1]
    WT_8[2.2] == t8[1.2]

    WT_9[2.1] == t9[1.1]
    WT_9[2.2] == t9[1.2]

    WT_10[2.1] == t10[1.1]
    WT_10[2.2] == t10[1.2]

    # Single WT
    # WT_1[2.1] == G_receiving[1.1] == Node1d
    # WT_1[2.2] == G_receiving[1.2] == Node1q
    # Two WTs
    # WT_1[2.1] == WT_2[2.1] == G_receiving[1.1] == Node1d
    # WT_1[2.2] == WT_2[2.2] == G_receiving[1.2] == Node1q
    # Ten WTs
    t1[2.1] == t2[2.1] == t3[2.1] ==t4[2.1] == G_receiving[1.1] == Node1d
    t1[2.2] == t2[2.2] == t3[2.2] ==t4[2.2] == G_receiving[1.2] == Node1q
end

omega_min = -1
omega_max = 3
omega_points = 1000

imp_WT_1, omega = determine_impedance(net, elim_elements = [:G_receiving], input_pins = Any[:Node1d, :Node1q], output_pins = Any[:gndd, :gndq], omega_range = (omega_min, omega_max, omega_points)) 

writedlm("./files/imp_WT_1.csv",  imp_WT_1, ',')
writedlm("./files/w_WT_1.csv",  omega, ',')


