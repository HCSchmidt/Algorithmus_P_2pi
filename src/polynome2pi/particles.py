from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class Particle:
    key: str
    name: str
    symbol: str
    theory_E: float
    sd_minus: float
    sd_plus: float
    theory_text: str | None = None

    @property
    def theory_str(self) -> str:
        return self.theory_text if self.theory_text else str(self.theory_E)


def get_particles() -> Dict[str, Particle]:
    particles = [
        Particle("electron", "e", "e", 1.0, -0.005, 0.0, "1.00000000000(31)"),
        Particle("u", "u", "u", 4.18, -0.51, 0.96, "4.18(-0.51)(0.96)"),
        Particle("d", "d", "d", 9.14, -0.33, 0.94, "9.14(-0.33)(0.94)"),
        Particle("s", "s", "s", 182.8, -6.6, 16.8, "182.8(-6.6)(16.8)"),
        Particle("muon", "Muon", "μ", 206.7682827, -0.0000046, 0.0000046, "206.7682827(46)"),
        Particle("pion_0", "Pion 0", "π0", 264.1430, -0.0009, 0.0009, "264.1430(9)"),
        Particle("pion_pm", "Pion +-", "π±", 273.13243, -0.00035, 0.00035, "273.13243(35)"),
        Particle("k_pm", "K +-", "K±", 966.102, -0.021, 0.021, "966.102(21)"),
        Particle("kl_0", "KL 0", "K_L0", 973.800, -0.026, 0.026, "973.800(26)"),
        Particle("ks_0", "KS 0", "K_S0", 973.800, -0.026, 0.026, "973.800(26)"),
        Particle("eta", "Eta", "η", 1072.139, -0.035, 0.035, "1072.139(35)"),
        Particle("rho_pm", "Rho +-", "ρ±", 1506.0, -1.0, 1.0, "1506(1)"),
        Particle("rho_0", "Rho 0", "ρ0", 1517.14, -0.49, 0.49, "1517.14(49)"),
        Particle("omega", "Omega", "ω", 1531.62, -0.25, 0.25, "1531.62(25)"),
        Particle("kstar_pm", "K* +-", "K*±", 1745.2, -0.1, 0.1, "1745.2(1)"),
        Particle("kstar_0", "K* 0", "K*0", 1752.6, -0.1, 0.1, "1752.6(1)"),
        Particle(
            "proton", "Proton", "p", 1836.152673426, -0.000000032, 0.000000032, "1836.152673426(32)"
        ),
        Particle("h_atom", "H", "H", 1837.47, -0.29, 0.20, "1837.47(-0.29)(0.20)"),
        Particle(
            "neutron", "Neutron", "n", 1838.68366200, -0.00000074, 0.00000074, "1838.68366200(74)"
        ),
        Particle("eta_prime", "Eta`", "η′", 1874.32, -0.11, 0.11, "1874.32(11)"),
        Particle("phi", "Phi", "φ", 1995.035, -0.031, 0.031, "1995.035(31)"),
        Particle("c", "c", "c", 2485.0, -39.0, 39.0, "2485(-39)(39)"),
        Particle("tau", "Tau", "τ", 3477.23, -0.23, 0.23, "3477.23(23)"),
    ]
    return {p.key: p for p in particles}
