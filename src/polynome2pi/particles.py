from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class Particle:
    """
    Physical particle definition.

    This object represents a physical constant in the P(2π) model.
    Instances must be immutable.
    """
    key: str                 # stable internal identifier (used in code)
    name: str                # display name (human-readable)
    symbol: str              # short label for plots
    theory_E: float          # theoretical energy (in m_e units)
    sd_minus: float          # lower deviation
    sd_plus: float           # upper deviation

    # Optional metadata (kept for reporting / future use)
    mass_text: str = ""
    charge: str = ""
    spin: str = ""
    parity: str = ""


# ---------------------------------------------------------------------
# Particle registry (single source of truth)
# ---------------------------------------------------------------------

PARTICLES: Dict[str, Particle] = {
    # -----------------------------------------------------------------
    # Leptons
    # -----------------------------------------------------------------
    "electron": Particle(
        key="electron",
        name="e",
        symbol="e",
        theory_E=1.0,
        sd_minus=-0.005,
        sd_plus=0.0,
        mass_text="1.00000000000(31)",
        charge="-1",
        spin="1/2",
    ),

    "muon": Particle(
        key="muon",
        name="Muon",
        symbol="μ",
        theory_E=206.7682827,
        sd_minus=-0.0000046,
        sd_plus=0.0000046,
        mass_text="206.7682827(46)",
        charge="-1",
        spin="1/2",
    ),

    "tau": Particle(
        key="tau",
        name="Tau",
        symbol="τ",
        theory_E=3477.23,
        sd_minus=-0.23,
        sd_plus=0.23,
        mass_text="3477.23(23)",
        charge="-1",
        spin="1/2",
    ),

    # -----------------------------------------------------------------
    # Quarks (effective masses used in the model)
    # -----------------------------------------------------------------
    "u": Particle(
        key="u",
        name="u",
        symbol="u",
        theory_E=4.18,
        sd_minus=-0.51,
        sd_plus=0.96,
        mass_text="4.18(-0.51)(0.96)",
        charge="+2/3",
    ),

    "d": Particle(
        key="d",
        name="d",
        symbol="d",
        theory_E=9.14,
        sd_minus=-0.33,
        sd_plus=0.94,
        mass_text="9.14(-0.33)(0.94)",
        charge="-1/3",
    ),

    "s": Particle(
        key="s",
        name="s",
        symbol="s",
        theory_E=182.8,
        sd_minus=-6.6,
        sd_plus=16.8,
        mass_text="182.8(-6.6)(16.8)",
        charge="-1/3",
    ),

    "c": Particle(
        key="c",
        name="c",
        symbol="c",
        theory_E=2485.0,
        sd_minus=-39.0,
        sd_plus=39.0,
        mass_text="2485(-39)(39)",
        charge="+2/3",
    ),

    # -----------------------------------------------------------------
    # Mesons
    # -----------------------------------------------------------------
    "pion_0": Particle(
        key="pion_0",
        name="Pion 0",
        symbol="π0",
        theory_E=264.1430,
        sd_minus=-0.0009,
        sd_plus=0.0009,
        mass_text="264.1430(9)",
    ),

    "pion_pm": Particle(
        key="pion_pm",
        name="Pion +-",
        symbol="π±",
        theory_E=273.13243,
        sd_minus=-0.00035,
        sd_plus=0.00035,
        mass_text="273.13243(35)",
    ),

    "rho_pm": Particle(
        key="rho_pm",
        name="Rho +-",
        symbol="ρ±",
        theory_E=1506.0,
        sd_minus=-1.0,
        sd_plus=1.0,
        mass_text="1506(1)",
    ),

    "rho_0": Particle(
        key="rho_0",
        name="Rho 0",
        symbol="ρ0",
        theory_E=1517.14,
        sd_minus=-0.49,
        sd_plus=0.49,
        mass_text="1517.14(49)",
    ),

    "omega": Particle(
        key="omega",
        name="Omega",
        symbol="ω",
        theory_E=1531.62,
        sd_minus=-0.25,
        sd_plus=0.25,
        mass_text="1531.62(25)",
    ),

    "eta": Particle(
        key="eta",
        name="Eta",
        symbol="η",
        theory_E=1072.139,
        sd_minus=-0.035,
        sd_plus=0.035,
        mass_text="1072.139(35)",
    ),

    "eta_prime": Particle(
        key="eta_prime",
        name="Eta`",
        symbol="η′",
        theory_E=1874.32,
        sd_minus=-0.11,
        sd_plus=0.11,
        mass_text="1874.32(11)",
    ),

    "phi": Particle(
        key="phi",
        name="Phi",
        symbol="φ",
        theory_E=1995.035,
        sd_minus=-0.031,
        sd_plus=0.031,
        mass_text="1995.035(31)",
    ),

    "kstar_pm": Particle(
        key="kstar_pm",
        name="K* +-",
        symbol="K*±",
        theory_E=1745.2,
        sd_minus=-0.1,
        sd_plus=0.1,
        mass_text="1745.2(1)",
    ),

    "kstar_0": Particle(
        key="kstar_0",
        name="K* 0",
        symbol="K*0",
        theory_E=1752.6,
        sd_minus=-0.1,
        sd_plus=0.1,
        mass_text="1752.6(1)",
    ),

    # -----------------------------------------------------------------
    # Baryons
    # -----------------------------------------------------------------
    "proton": Particle(
        key="proton",
        name="Proton",
        symbol="p",
        theory_E=1836.152673426,
        sd_minus=-0.000000032,
        sd_plus=0.000000032,
        mass_text="1836.152673426(32)",
        charge="+1",
        spin="1/2",
    ),

    "neutron": Particle(
        key="neutron",
        name="Neutron",
        symbol="n",
        theory_E=1838.68366200,
        sd_minus=-0.00000074,
        sd_plus=0.00000074,
        mass_text="1838.68366200(74)",
        spin="1/2",
    ),

    "deuteron": Particle(
        key="deuteron",
        name="Deuteron",
        symbol="D",
        theory_E=3670.4829677,
        sd_minus=-0.0000011,
        sd_plus=0.0000011,
        mass_text="3670.4829677(11)",
    ),

    # -----------------------------------------------------------------
    # Higgs / heavy
    # -----------------------------------------------------------------
    "higgs": Particle(
        key="higgs",
        name="Higgs",
        symbol="H",
        theory_E=244830.0,
        sd_minus=-210.0,
        sd_plus=210.0,
        mass_text="244830(210)",
    ),
}