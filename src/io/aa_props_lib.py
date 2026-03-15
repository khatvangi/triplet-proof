# src/io/aa_props_lib.py
# 24-ish properties for the 20 standard amino acids (1-letter codes).
# Values are canonical/approx; downstream we z-score + covariance-correct.

AA = {}

def add(aa, **kvs):
    AA[aa] = kvs

# Kyte–Doolittle hydropathy ~ hydro_kd
# Zamyatnin volume (Å^3) ~ vol_zam
# Grantham polarity ~ polarity_gr
# pI ~ isoelectric point
# sc_pka ~ side-chain pKa (0 if none)
# charge_pH7 ~ net side-chain charge at pH 7 (His ~ +0.1)
# donors/acceptors ~ side-chain HBD/HBA counts
# arom (0/1), ring_count, rotb (rotatable bonds)
# refractivity, bulkiness (Zimmerman)
# Chou–Fasman: cf_alpha, cf_beta, cf_turn
# acc_pref ~ relative solvent accessibility
# flex ~ normalized flexibility index
# helical_moment, beta_moment ~ simple moments (proxy)
# hphob_fr ~ hydrophobic surface tendency (proxy)
# sc_mass ~ side-chain mass (Da, approx minus backbone)

add("A", hydro_kd=1.8, vol_zam=88.6, polarity_gr=8.1, pI=6.0, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=0, rotb=0, refractivity=1.41, bulkiness=11.5,
    cf_alpha=1.42, cf_beta=0.83, cf_turn=0.66, acc_pref=0.35, flex=0.45,
    helical_moment=0.18, beta_moment=0.10, hphob_fr=0.62, sc_mass=15.0)

add("R", hydro_kd=-4.5, vol_zam=173.4, polarity_gr=10.5, pI=10.8, sc_pka=12.5, charge_pH7=+1,
    donors=4, acceptors=1, arom=0, ring_count=0, rotb=4, refractivity=3.04, bulkiness=14.3,
    cf_alpha=0.98, cf_beta=0.93, cf_turn=0.95, acc_pref=0.75, flex=0.54,
    helical_moment=0.05, beta_moment=0.08, hphob_fr=0.10, sc_mass=101.0)

add("N", hydro_kd=-3.5, vol_zam=114.1, polarity_gr=11.6, pI=5.4, sc_pka=0.0, charge_pH7=0,
    donors=1, acceptors=1, arom=0, ring_count=0, rotb=2, refractivity=1.93, bulkiness=12.3,
    cf_alpha=0.67, cf_beta=0.89, cf_turn=1.56, acc_pref=0.70, flex=0.51,
    helical_moment=0.07, beta_moment=0.12, hphob_fr=0.18, sc_mass=58.0)

add("D", hydro_kd=-3.5, vol_zam=111.1, polarity_gr=13.0, pI=2.8, sc_pka=3.9, charge_pH7=-1,
    donors=0, acceptors=1, arom=0, ring_count=0, rotb=2, refractivity=1.87, bulkiness=13.0,
    cf_alpha=1.01, cf_beta=0.54, cf_turn=1.46, acc_pref=0.80, flex=0.50,
    helical_moment=0.08, beta_moment=0.09, hphob_fr=0.12, sc_mass=59.0)

add("C", hydro_kd=2.5, vol_zam=108.5, polarity_gr=5.5, pI=5.1, sc_pka=8.3, charge_pH7=0,
    donors=0, acceptors=1, arom=0, ring_count=0, rotb=1, refractivity=1.76, bulkiness=13.5,
    cf_alpha=0.70, cf_beta=1.19, cf_turn=1.19, acc_pref=0.45, flex=0.46,
    helical_moment=0.16, beta_moment=0.12, hphob_fr=0.66, sc_mass=47.0)

add("Q", hydro_kd=-3.5, vol_zam=143.8, polarity_gr=10.5, pI=5.7, sc_pka=0.0, charge_pH7=0,
    donors=1, acceptors=1, arom=0, ring_count=0, rotb=3, refractivity=2.28, bulkiness=14.5,
    cf_alpha=1.11, cf_beta=1.10, cf_turn=0.98, acc_pref=0.68, flex=0.53,
    helical_moment=0.06, beta_moment=0.11, hphob_fr=0.16, sc_mass=72.0)

add("E", hydro_kd=-3.5, vol_zam=138.4, polarity_gr=12.3, pI=3.2, sc_pka=4.2, charge_pH7=-1,
    donors=0, acceptors=1, arom=0, ring_count=0, rotb=3, refractivity=2.10, bulkiness=13.6,
    cf_alpha=1.51, cf_beta=0.37, cf_turn=0.74, acc_pref=0.85, flex=0.52,
    helical_moment=0.04, beta_moment=0.07, hphob_fr=0.13, sc_mass=73.0)

add("G", hydro_kd=-0.4, vol_zam=60.1, polarity_gr=9.0, pI=6.0, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=0, rotb=0, refractivity=0.90, bulkiness=3.4,
    cf_alpha=0.57, cf_beta=0.75, cf_turn=1.64, acc_pref=0.60, flex=0.70,
    helical_moment=0.12, beta_moment=0.12, hphob_fr=0.30, sc_mass=1.0)

add("H", hydro_kd=-3.2, vol_zam=153.2, polarity_gr=10.4, pI=7.6, sc_pka=6.0, charge_pH7=+0.1,
    donors=1, acceptors=1, arom=1, ring_count=1, rotb=2, refractivity=2.55, bulkiness=13.7,
    cf_alpha=1.00, cf_beta=0.87, cf_turn=0.95, acc_pref=0.65, flex=0.49,
    helical_moment=0.07, beta_moment=0.10, hphob_fr=0.20, sc_mass=82.0)

add("I", hydro_kd=4.5, vol_zam=166.7, polarity_gr=5.2, pI=6.0, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=0, rotb=2, refractivity=2.66, bulkiness=19.8,
    cf_alpha=1.08, cf_beta=1.60, cf_turn=0.47, acc_pref=0.25, flex=0.45,
    helical_moment=0.24, beta_moment=0.06, hphob_fr=0.85, sc_mass=57.0)

add("L", hydro_kd=3.8, vol_zam=166.7, polarity_gr=4.9, pI=6.0, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=0, rotb=2, refractivity=2.59, bulkiness=21.4,
    cf_alpha=1.21, cf_beta=1.30, cf_turn=0.59, acc_pref=0.28, flex=0.45,
    helical_moment=0.23, beta_moment=0.07, hphob_fr=0.83, sc_mass=57.0)

add("K", hydro_kd=-3.9, vol_zam=168.6, polarity_gr=11.3, pI=9.7, sc_pka=10.5, charge_pH7=+1,
    donors=2, acceptors=1, arom=0, ring_count=0, rotb=4, refractivity=2.93, bulkiness=15.7,
    cf_alpha=1.16, cf_beta=0.74, cf_turn=1.01, acc_pref=0.73, flex=0.53,
    helical_moment=0.05, beta_moment=0.08, hphob_fr=0.12, sc_mass=72.0)

add("M", hydro_kd=1.9, vol_zam=162.9, polarity_gr=5.7, pI=5.7, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=1, arom=0, ring_count=0, rotb=3, refractivity=2.47, bulkiness=16.3,
    cf_alpha=1.45, cf_beta=1.05, cf_turn=0.60, acc_pref=0.40, flex=0.48,
    helical_moment=0.18, beta_moment=0.10, hphob_fr=0.70, sc_mass=75.0)

add("F", hydro_kd=2.8, vol_zam=189.9, polarity_gr=5.2, pI=5.5, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=1, ring_count=1, rotb=2, refractivity=3.02, bulkiness=19.8,
    cf_alpha=1.13, cf_beta=1.38, cf_turn=0.60, acc_pref=0.30, flex=0.44,
    helical_moment=0.20, beta_moment=0.08, hphob_fr=0.88, sc_mass=91.0)

add("P", hydro_kd=-1.6, vol_zam=112.7, polarity_gr=8.0, pI=6.3, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=1, rotb=0, refractivity=1.77, bulkiness=17.4,
    cf_alpha=0.57, cf_beta=0.55, cf_turn=1.52, acc_pref=0.55, flex=0.36,
    helical_moment=0.08, beta_moment=0.15, hphob_fr=0.35, sc_mass=42.0)

add("S", hydro_kd=-0.8, vol_zam=89.0, polarity_gr=9.2, pI=5.7, sc_pka=13.0, charge_pH7=0,
    donors=1, acceptors=1, arom=0, ring_count=0, rotb=1, refractivity=1.43, bulkiness=9.2,
    cf_alpha=0.77, cf_beta=0.75, cf_turn=1.43, acc_pref=0.75, flex=0.51,
    helical_moment=0.11, beta_moment=0.12, hphob_fr=0.28, sc_mass=31.0)

add("T", hydro_kd=-0.7, vol_zam=116.1, polarity_gr=8.6, pI=5.6, sc_pka=13.0, charge_pH7=0,
    donors=1, acceptors=1, arom=0, ring_count=0, rotb=2, refractivity=1.76, bulkiness=15.7,
    cf_alpha=0.83, cf_beta=1.19, cf_turn=1.21, acc_pref=0.70, flex=0.49,
    helical_moment=0.12, beta_moment=0.11, hphob_fr=0.32, sc_mass=45.0)

add("W", hydro_kd=-0.9, vol_zam=227.8, polarity_gr=5.4, pI=5.9, sc_pka=0.0, charge_pH7=0,
    donors=1, acceptors=1, arom=1, ring_count=2, rotb=2, refractivity=3.75, bulkiness=21.6,
    cf_alpha=1.08, cf_beta=1.37, cf_turn=0.57, acc_pref=0.25, flex=0.42,
    helical_moment=0.22, beta_moment=0.07, hphob_fr=0.86, sc_mass=130.0)

add("Y", hydro_kd=-1.3, vol_zam=193.6, polarity_gr=6.2, pI=5.7, sc_pka=10.1, charge_pH7=0,
    donors=1, acceptors=1, arom=1, ring_count=1, rotb=2, refractivity=3.15, bulkiness=18.0,
    cf_alpha=0.69, cf_beta=1.47, cf_turn=1.14, acc_pref=0.35, flex=0.45,
    helical_moment=0.19, beta_moment=0.09, hphob_fr=0.80, sc_mass=107.0)

add("V", hydro_kd=4.2, vol_zam=140.0, polarity_gr=5.9, pI=6.0, sc_pka=0.0, charge_pH7=0,
    donors=0, acceptors=0, arom=0, ring_count=0, rotb=1, refractivity=2.35, bulkiness=21.6,
    cf_alpha=1.06, cf_beta=1.70, cf_turn=0.50, acc_pref=0.27, flex=0.46,
    helical_moment=0.24, beta_moment=0.06, hphob_fr=0.84, sc_mass=43.0)

# handy ordered list of column names
def aa_property_columns():
    return [
        "hydro_kd","vol_zam","polarity_gr","pI","sc_pka","charge_pH7",
        "donors","acceptors","arom","ring_count","rotb","refractivity",
        "bulkiness","cf_alpha","cf_beta","cf_turn","acc_pref","flex",
        "helical_moment","beta_moment","hphob_fr","sc_mass"
    ]

