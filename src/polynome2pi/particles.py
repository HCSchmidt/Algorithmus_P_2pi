# src/polynome2pi/particles.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class Particle:
    """
    Particle definition used for:
      - matching scan results against theory value + tolerances (sd_minus/sd_plus)
      - writing legacy CSV (theory_text keeps your old strings like "4.18(-0.51)(0.96)")

    Conventions:
      - theory_E is numeric center value used for comparisons
      - sd_minus / sd_plus are numeric deltas (can be negative/positive)
      - theory_text is the exact display string from the legacy table (optional)
    """
    key: str
    name: str
    symbol: str
    theory_E: float
    sd_minus: float = 0.0
    sd_plus: float = 0.0
    theory_text: str | None = None

    @property
    def theory_str(self) -> str:
        return self.theory_text if self.theory_text else str(self.theory_E)


# IMPORTANT:
# - Keys should stay stable (engine/results dict keying).
# - Order matters for stable plotting colors if you map by insertion order.
PARTICLES: Dict[str, Particle] = {
    # leptons / quarks
    "electron": Particle(
        key="electron",
        name="e",
        symbol="e",
        theory_E=1.0,
        sd_minus=-0.005,
        sd_plus=0.000,
        theory_text="1.00000000000(31)",
    ),
    "u": Particle(
        key="u",
        name="u",
        symbol="u",
        theory_E=4.18,
        sd_minus=-0.51,
        sd_plus=0.96,
        theory_text="4.18(-0.51)(0.96)",
    ),
    "d": Particle(
        key="d",
        name="d",
        symbol="d",
        theory_E=9.14,
        sd_minus=-0.33,
        sd_plus=0.94,
        theory_text="9.14(-0.33)(0.94)",
    ),
    "s": Particle(
        key="s",
        name="s",
        symbol="s",
        theory_E=182.8,
        sd_minus=-6.6,
        sd_plus=16.8,
        theory_text="182.8(-6.6)(16.8)",
    ),
    "muon": Particle(
        key="muon",
        name="Muon",
        symbol="μ",
        theory_E=206.7682827,
        sd_minus=-0.0000046,
        sd_plus=0.0000046,
        theory_text="206.7682827(46)",
    ),

    # mesons
    "pion_0": Particle(
        key="pion_0",
        name="Pion 0",
        symbol="π0",
        theory_E=264.1430,
        sd_minus=-0.0009,
        sd_plus=0.0009,
        theory_text="264.1430(9)",
    ),
    "pion_pm": Particle(
        key="pion_pm",
        name="Pion +-",
        symbol="π±",
        theory_E=273.13243,
        sd_minus=-0.00035,
        sd_plus=0.00035,
        theory_text="273.13243(35)",
    ),
    "k_pm": Particle(
        key="k_pm",
        name="K +-",
        symbol="K±",
        theory_E=966.102,
        sd_minus=-0.021,
        sd_plus=0.021,
        theory_text="966.102(21)",
    ),
    "kl_0": Particle(
        key="kl_0",
        name="KL 0",
        symbol="K_L0",
        theory_E=973.800,
        sd_minus=-0.026,
        sd_plus=0.026,
        theory_text="973.800(26)",
    ),
    "ks_0": Particle(
        key="ks_0",
        name="KS 0",
        symbol="K_S0",
        theory_E=973.800,
        sd_minus=-0.026,
        sd_plus=0.026,
        theory_text="973.800(26)",
    ),
    "eta": Particle(
        key="eta",
        name="Eta",
        symbol="η",
        theory_E=1072.139,
        sd_minus=-0.035,
        sd_plus=0.035,
        theory_text="1072.139(35)",
    ),
    "rho_pm": Particle(
        key="rho_pm",
        name="Rho +-",
        symbol="ρ±",
        theory_E=1506.0,
        sd_minus=-1.0,
        sd_plus=1.0,
        theory_text="1506(1)",
    ),
    "rho_0": Particle(
        key="rho_0",
        name="Rho 0",
        symbol="ρ0",
        theory_E=1517.14,
        sd_minus=-0.49,
        sd_plus=0.49,
        theory_text="1517.14(49)",
    ),
    "omega": Particle(
        key="omega",
        name="Omega",
        symbol="ω",
        theory_E=1531.62,
        sd_minus=-0.25,
        sd_plus=0.25,
        theory_text="1531.62(25)",
    ),
    "kstar_pm": Particle(
        key="kstar_pm",
        name="K* +-",
        symbol="K*±",
        theory_E=1745.2,
        sd_minus=-0.1,
        sd_plus=0.1,
        theory_text="1745.2(1)",
    ),
    "kstar_0": Particle(
        key="kstar_0",
        name="K* 0",
        symbol="K*0",
        theory_E=1752.6,
        sd_minus=-0.1,
        sd_plus=0.1,
        theory_text="1752.6(1)",
    ),

    # baryons / atoms
    "proton": Particle(
        key="proton",
        name="Proton",
        symbol="p",
        theory_E=1836.152673426,
        sd_minus=-0.000000032,
        sd_plus=0.000000032,
        theory_text="1836.152673426(32)",
    ),
    "h_atom": Particle(
        key="h_atom",
        name="H",
        symbol="H",
        theory_E=1837.47,
        sd_minus=-0.29,
        sd_plus=0.20,
        theory_text="1837.47(-0.29)(0.20)",
    ),
    "neutron": Particle(
        key="neutron",
        name="Neutron",
        symbol="n",
        theory_E=1838.68366200,
        sd_minus=-0.00000074,
        sd_plus=0.00000074,
        theory_text="1838.68366200(74)",
    ),

    # more mesons
    "eta_prime": Particle(
        key="eta_prime",
        name="Eta`",
        symbol="η′",
        theory_E=1874.32,
        sd_minus=-0.11,
        sd_plus=0.11,
        theory_text="1874.32(11)",
    ),
    "phi": Particle(
        key="phi",
        name="Phi",
        symbol="φ",
        theory_E=1995.035,
        sd_minus=-0.031,
        sd_plus=0.031,
        theory_text="1995.035(31)",
    ),

    # heavy flavors / tau
    "c": Particle(
        key="c",
        name="c",
        symbol="c",
        theory_E=2485.0,
        sd_minus=-39.0,
        sd_plus=39.0,
        theory_text="2485(-39)(39)",
    ),
    "tau": Particle(
        key="tau",
        name="Tau",
        symbol="τ",
        theory_E=3477.23,
        sd_minus=-0.23,
        sd_plus=0.23,
        theory_text="3477.23(23)",
    ),

    # D mesons / deuteron
    "d0": Particle(
        key="d0",
        name="D 0",
        symbol="D0",
        theory_E=3649.38,
        sd_minus=-0.10,
        sd_plus=0.10,
        theory_text="3649.38(10)",
    ),
    "dplus": Particle(
        key="dplus",
        name="D +",
        symbol="D+",
        theory_E=3658.81,
        sd_minus=-0.10,
        sd_plus=0.10,
        theory_text="3658.81(10)",
    ),
    "deuteron": Particle(
        key="deuteron",
        name="Deuteron",
        symbol="D",
        theory_E=3670.4829677,
        sd_minus=-0.0000011,
        sd_plus=0.0000011,
        theory_text="3670.4829677(11)",
    ),
    "ds_plus": Particle(
        key="ds_plus",
        name="DS +",
        symbol="D_s+",
        theory_E=3851.94,
        sd_minus=-0.13,
        sd_plus=0.13,
        theory_text="3851.94(13)",
    ),

    # Higgs / top (kept for completeness)
    "higgs": Particle(
        key="higgs",
        name="Higgs",
        symbol="H",
        theory_E=244830.0,
        sd_minus=-210.0,
        sd_plus=210.0,
        theory_text="244830(210)",
    ),
    "t": Particle(
        key="t",
        name="t",
        symbol="t",
        theory_E=337710.0,
        sd_minus=-570.0,
        sd_plus=570.0,
        theory_text="337710(570)",
    ),
}