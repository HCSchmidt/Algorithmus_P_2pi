from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional


# -----------------------------------------------------------------------------
# Physics particles
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Particle:
    """One entry used by the scan.

    The scan only needs:
      - key, name, symbol
      - theory_E (in units of m_e)
      - sd_minus / sd_plus (tolerance window around theory_E)

    Everything else is optional metadata (for labels, exports, etc.).
    """

    # Required by the scan engine
    key: str
    name: str
    symbol: str
    theory_E: float
    sd_minus: float
    sd_plus: float

    # Optional display / metadata
    category: str = ""
    theory_text: Optional[str] = None
    half_life_s: str = ""
    charge: str = ""
    spin: str = ""
    parity: str = ""
    composition: str = ""  # e.g. quark content (string/LaTeX)

    # Controls whether this entry participates in the scan by default.
    # (Massless bosons would otherwise match lots of low-energy points.)
    include_in_scan: bool = True

    @property
    def theory_str(self) -> str:
        return self.theory_text if self.theory_text else str(self.theory_E)


def _particle_rows() -> List[Particle]:
    """Single source of truth for particle definitions."""

    return [
        # --- massless bosons / gluons (kept for reference; excluded from scan by default)
        Particle(
            key="photon",
            name="Photon",
            symbol="γ",
            theory_E=0.0,
            sd_minus=0.0,
            sd_plus=0.0,
            category="Boson",
            theory_text="0",
            half_life_s="∞",
            charge="0",
            spin="1",
            parity="-",
            composition="",
            include_in_scan=False,
        ),
        Particle(
            key="gluon_r",
            name="r",
            symbol="r",
            theory_E=0.0,
            sd_minus=0.0,
            sd_plus=0.0,
            category="Gluon",
            theory_text="0",
            half_life_s="∞",
            charge="0",
            spin="1",
            parity="-",
            composition="",
            include_in_scan=False,
        ),
        Particle(
            key="gluon_b",
            name="b",
            symbol="b",
            theory_E=0.0,
            sd_minus=0.0,
            sd_plus=0.0,
            category="Gluon",
            theory_text="0",
            half_life_s="∞",
            charge="0",
            spin="1",
            parity="-",
            composition="",
            include_in_scan=False,
        ),
        Particle(
            key="gluon_g",
            name="g",
            symbol="g",
            theory_E=0.0,
            sd_minus=0.0,
            sd_plus=0.0,
            category="Gluon",
            theory_text="0",
            half_life_s="∞",
            charge="0",
            spin="1",
            parity="-",
            composition="",
            include_in_scan=False,
        ),

        # --- leptons / quarks
        Particle(
            key="electron",
            name="e",
            symbol="e",
            theory_E=1.0,
            sd_minus=-0.005,
            sd_plus=0.0,
            category="Lepton",
            theory_text="1.00000000000(31)",
            half_life_s="",
            charge="-1",
            spin="1/2",
            parity="",
            composition="e",
        ),
        Particle(
            key="u",
            name="u",
            symbol="u",
            theory_E=4.18,
            sd_minus=-0.51,
            sd_plus=0.96,
            category="Quark",
            theory_text="4.18(-0.51)(0.96)",
            half_life_s="",
            charge="+2/3",
            spin="1/2",
            parity="",
            composition="u",
        ),
        Particle(
            key="d",
            name="d",
            symbol="d",
            theory_E=9.14,
            sd_minus=-0.33,
            sd_plus=0.94,
            category="Quark",
            theory_text="9.14(-0.33)(0.94)",
            half_life_s="",
            charge="-1/3",
            spin="1/2",
            parity="",
            composition="d",
        ),
        Particle(
            key="s",
            name="s",
            symbol="s",
            theory_E=182.8,
            sd_minus=-6.6,
            sd_plus=16.8,
            category="Quark",
            theory_text="182.8(-6.6)(16.8)",
            half_life_s="",
            charge="-1/3",
            spin="1/2",
            parity="",
            composition="s",
        ),
        Particle(
            key="muon",
            name="Muon",
            symbol="μ",
            theory_E=206.7682827,
            sd_minus=-0.0000046,
            sd_plus=0.0000046,
            category="Lepton",
            theory_text="206.7682827(46)",
            half_life_s="2.1969811(22)e-6",
            charge="-1",
            spin="1/2",
            parity="",
            composition="muon",
        ),

        # --- mesons
        Particle(
            key="pion_0",
            name="Pion 0",
            symbol="π0",
            theory_E=264.1430,
            sd_minus=-0.0009,
            sd_plus=0.0009,
            category="Meson",
            theory_text="264.1430(9)",
            half_life_s="8.52(18)e-17",
            charge="0",
            spin="0",
            parity="-",
            composition=r"$u\overline{d}-\overline{u}d$",
        ),
        Particle(
            key="pion_pm",
            name="Pion +-",
            symbol="π±",
            theory_E=273.13243,
            sd_minus=-0.00035,
            sd_plus=0.00035,
            category="Meson",
            theory_text="273.13243(35)",
            half_life_s="2.6033(5)e-8",
            charge="±1",
            spin="0",
            parity="-",
            composition=r"$u\overline{u},\overline{d}d$",
        ),
        Particle(
            key="k_pm",
            name="K +-",
            symbol="K±",
            theory_E=966.102,
            sd_minus=-0.021,
            sd_plus=0.021,
            category="Meson",
            theory_text="966.102(21)",
            half_life_s="1.2380(20)e-8",
            charge="±1",
            spin="0",
            parity="-",
            composition=r"$u\overline{s},s\overline{u}$",
        ),
        Particle(
            key="kl_0",
            name="KL 0",
            symbol="K_L0",
            theory_E=973.800,
            sd_minus=-0.026,
            sd_plus=0.026,
            category="Meson",
            theory_text="973.800(26)",
            half_life_s="5.116(21)e-8",
            charge="0",
            spin="0",
            parity="-",
            composition=r"$d\overline{s},s\overline{d}$",
        ),
        Particle(
            key="ks_0",
            name="KS 0",
            symbol="K_S0",
            theory_E=973.800,
            sd_minus=-0.026,
            sd_plus=0.026,
            category="Meson",
            theory_text="973.800(26)",
            half_life_s="8.954(4)e-11",
            charge="0",
            spin="0",
            parity="",
            composition=r"$d\overline{s},s\overline{d}$",
        ),
        Particle(
            key="eta",
            name="Eta",
            symbol="η",
            theory_E=1072.139,
            sd_minus=-0.035,
            sd_plus=0.035,
            category="Meson",
            theory_text="1072.139(35)",
            half_life_s="5e-19",
            charge="0",
            spin="1/2",
            parity="-",
            composition=r"$u\overline{u}+\overline{d}d-2s\overline{s}$",
        ),
        Particle(
            key="rho_pm",
            name="Rho +-",
            symbol="ρ±",
            theory_E=1506.0,
            sd_minus=-1.0,
            sd_plus=1.0,
            category="Meson",
            theory_text="1506(1)",
            half_life_s="4e-24",
            charge="±1",
            spin="1",
            parity="-",
            composition=r"$u\overline{u},\overline{d}d$",
        ),
        Particle(
            key="rho_0",
            name="Rho 0",
            symbol="ρ0",
            theory_E=1517.14,
            sd_minus=-0.49,
            sd_plus=0.49,
            category="Meson",
            theory_text="1517.14(49)",
            half_life_s="4e-24",
            charge="0",
            spin="1",
            parity="-",
            composition=r"$u\overline{u}-\overline{d}d$",
        ),
        Particle(
            key="omega",
            name="Omega",
            symbol="ω",
            theory_E=1531.62,
            sd_minus=-0.25,
            sd_plus=0.25,
            category="Meson",
            theory_text="1531.62(25)",
            half_life_s="7.75(7)e-23",
            charge="0",
            spin="1",
            parity="-",
            composition=r"$u\overline{u}+\overline{d}d$",
        ),
        Particle(
            key="kstar_pm",
            name="K* +-",
            symbol="K*±",
            theory_E=1745.2,
            sd_minus=-0.1,
            sd_plus=0.1,
            category="Meson",
            theory_text="1745.2(1)",
            half_life_s="1.3e-23",
            charge="±1",
            spin="0",
            parity="",
            composition=r"$d\overline{s},s\overline{d}$",
        ),
        Particle(
            key="kstar_0",
            name="K* 0",
            symbol="K*0",
            theory_E=1752.6,
            sd_minus=-0.1,
            sd_plus=0.1,
            category="Meson",
            theory_text="1752.6(1)",
            half_life_s="1.3e-23",
            charge="0",
            spin="0",
            parity="",
            composition=r"$d\overline{s},s\overline{d}$",
        ),

        # --- nucleons / atoms
        Particle(
            key="proton",
            name="Proton",
            symbol="p",
            theory_E=1836.152673426,
            sd_minus=-0.000000032,
            sd_plus=0.000000032,
            category="Nukleon",
            theory_text="1836.152673426(32)",
            half_life_s="",
            charge="1",
            spin="1/2",
            parity="1",
            composition="uud",
        ),
        Particle(
            key="h_atom",
            name="H",
            symbol="H",
            theory_E=1837.47,
            sd_minus=-0.29,
            sd_plus=0.20,
            category="Atom",
            theory_text="1837.47(-0.29)(0.20)",
            half_life_s="",
            charge="0",
            spin="0",
            parity="",
            composition="H",
        ),
        Particle(
            key="neutron",
            name="Neutron",
            symbol="n",
            theory_E=1838.68366200,
            sd_minus=-0.00000074,
            sd_plus=0.00000074,
            category="Nukleon",
            theory_text="1838.68366200(74)",
            half_life_s="878.4(5)",
            charge="0",
            spin="1/2",
            parity="1",
            composition="udd",
        ),

        # --- heavier
        Particle(
            key="eta_prime",
            name="Eta`",
            symbol="η′",
            theory_E=1874.32,
            sd_minus=-0.11,
            sd_plus=0.11,
            category="Meson",
            theory_text="1874.32(11)",
            half_life_s="3.32(15)e-21",
            charge="0",
            spin="1/2",
            parity="-",
            composition=r"$u\overline{u}+\overline{d}d+s\overline{s}$",
        ),
        Particle(
            key="phi",
            name="Phi",
            symbol="φ",
            theory_E=1995.035,
            sd_minus=-0.031,
            sd_plus=0.031,
            category="Meson",
            theory_text="1995.035(31)",
            half_life_s="1.55(0.01)e-22",
            charge="0",
            spin="1",
            parity="-",
            composition=r"$s\overline{s}(most)$",
        ),
        Particle(
            key="c",
            name="c",
            symbol="c",
            theory_E=2485.0,
            sd_minus=-39.0,
            sd_plus=39.0,
            category="Quark",
            theory_text="2485(-39)(39)",
            half_life_s="",
            charge="+2/3",
            spin="1/2",
            parity="",
            composition="c",
        ),
        Particle(
            key="tau",
            name="Tau",
            symbol="τ",
            theory_E=3477.23,
            sd_minus=-0.23,
            sd_plus=0.23,
            category="Lepton",
            theory_text="3477.23(23)",
            half_life_s="290.3(5)e-15",
            charge="-1",
            spin="1/2",
            parity="",
            composition="tau",
        ),
        Particle(
            key="d_0",
            name="D 0",
            symbol="D0",
            theory_E=3649.38,
            sd_minus=-0.10,
            sd_plus=0.10,
            category="Meson",
            theory_text="3649.38(10)",
            half_life_s="4.101(15)e-13",
            charge="±1",
            spin="0",
            parity="-",
            composition=r"$c\overline{u},u\overline{c}$",
        ),
        Particle(
            key="d_plus",
            name="D +",
            symbol="D+",
            theory_E=3658.81,
            sd_minus=-0.10,
            sd_plus=0.10,
            category="Meson",
            theory_text="3658.81(10)",
            half_life_s="1.040(7)e-12",
            charge="±1",
            spin="0",
            parity="-",
            composition=r"$c\overline{d},d\overline{c}$",
        ),
        Particle(
            key="deuteron",
            name="Deuteron",
            symbol="D",
            theory_E=3670.4829677,
            sd_minus=-0.0000011,
            sd_plus=0.0000011,
            category="Nukleon",
            theory_text="3670.4829677(11)",
            half_life_s="",
            charge="0",
            spin="0",
            parity="",
            composition="D",
        ),
        Particle(
            key="ds_plus",
            name="DS +",
            symbol="D_s+",
            theory_E=3851.94,
            sd_minus=-0.13,
            sd_plus=0.13,
            category="Meson",
            theory_text="3851.94(13)",
            half_life_s="5.04(4)e-13",
            charge="±1",
            spin="0",
            parity="-",
            composition=r"$c\overline{s},s\overline{c}$",
        ),
        Particle(
            key="higgs",
            name="Higgs",
            symbol="Higgs",
            theory_E=244830.0,
            sd_minus=-210.0,
            sd_plus=210.0,
            category="Boson",
            theory_text="244830(210)",
            half_life_s="",
            charge="0",
            spin="0",
            parity="",
            composition="",
        ),
        Particle(
            key="t",
            name="t",
            symbol="t",
            theory_E=337710.0,
            sd_minus=-570.0,
            sd_plus=570.0,
            category="Quark",
            theory_text="337710(570)",
            half_life_s="",
            charge="+2/3",
            spin="1/2",
            parity="",
            composition="t",
        ),
    ]


def get_particles(*, include_reference_massless: bool = False) -> Dict[str, Particle]:
    """Return particles used by the scan.

    By default this returns the *scan* set (exclude massless reference bosons/gluons).
    If you really want them in the scan, pass include_reference_massless=True.
    """

    rows = _particle_rows()
    if include_reference_massless:
        selected = rows
    else:
        selected = [p for p in rows if p.include_in_scan]
    return {p.key: p for p in selected}


def iter_all_particles() -> Iterable[Particle]:
    """Iterate over *all* particle definitions (including excluded reference rows)."""

    yield from _particle_rows()


# -----------------------------------------------------------------------------
# Solar system bodies (kept here because they came from the same legacy tables)
# -----------------------------------------------------------------------------


AE = 0.387098  # Astronomical unit scale used in the legacy table (as given)


@dataclass(frozen=True)
class SolarBody:
    category: str
    name: str
    radius_km: float
    polar_radius_km: float
    tilt_deg: float
    orbit_radius_km: float
    orbit_radius_minus_sd: float
    orbit_radius_plus_sd: float
    period_days: float
    ecliptic_deg: float


def get_solar_system_bodies() -> Dict[str, SolarBody]:
    """Returns the "Sun table" entries as structured data.

    Note: values are kept as-is from the provided legacy snippet.
    """

    bodies: List[SolarBody] = [
        SolarBody("Sun", "Sun", 696_342, 0.0, 0.0, 0.0, 0.0, 0.0, 25.38, 0.0),
        SolarBody(
            "Planet",
            "Mercury",
            4881 / 2,
            4876.6 / 2,
            0.034,
            57_900_000,
            0.308 * AE,
            0.467 * AE,
            87.969,
            7.004,
        ),
        SolarBody(
            "Planet",
            "Venus",
            12103.6 / 2,
            12103.6 / 2,
            177.36,
            108_200_000,
            0.718 * AE,
            0.728 * AE,
            224.701,
            3.395,
        ),
        SolarBody(
            "Planet",
            "Earth",
            6_378_137.0,
            6_356_752.314,
            23.44,
            149_600_000,
            0.983 * AE,
            1.017 * AE,
            365.256,
            0.0001,
        ),
        SolarBody(
            "Planet",
            "Mars",
            6792.4 / 2,
            6752.4 / 2,
            25.19,
            227_990_000,
            1.382 * AE,
            1.666 * AE,
            686.980,
            1.8506,
        ),
        SolarBody(
            "Moon",
            "Moon",
            3474 / 2,
            3474 / 2,
            6.68,
            384_400,
            363_300,
            405_500,
            27.3217,
            5.145,
        ),
    ]

    # key by lowercase name
    return {b.name.lower(): b for b in bodies}
