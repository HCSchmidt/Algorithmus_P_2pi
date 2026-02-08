from polynome2pi.energy_model import EnergyModel
from polynome2pi.engine import ScanEngine
from polynome2pi.particles import get_particles
from polynome2pi.presets import preset_for_sector
from polynome2pi.scan.plotting import plot_scan
from polynome2pi.scan.report import write_results_csv


def run_scan(sector, results_dir):
    model = EnergyModel()
    preset = preset_for_sector(sector)

    engine = ScanEngine(preset=preset, model=model)
    particles = get_particles()
    outputs = engine.run(particles)

    base_name = f"scan_{sector.value}"
    png_path = results_dir / f"{base_name}.png"
    csv_path = results_dir / f"{base_name}.csv"
    title = (f"P(2π) scan – sector: {preset.name}",)
    plot_scan(
        out_png=png_path,
        particles=particles,
        matched_points=outputs.matched_points,
        unmatched_segments=outputs.unmatched_segments,
        title=title,
    )

    write_results_csv(
        path=csv_path,
        sector=sector,
        particles=particles,
        bins_by_particle=outputs.bins_by_particle,
    )

    return png_path
