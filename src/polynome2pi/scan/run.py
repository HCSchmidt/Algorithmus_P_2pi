from polynome2pi.energy_model import EnergyModel
from polynome2pi.engine import ScanEngine
from polynome2pi.particles import get_particles
from polynome2pi.presets import preset_for_sector
from polynome2pi.scan.plotting import plot_scan, plot_match_grid_3d_scatter
from polynome2pi.scan.report import write_results_csv, write_results_txt


def run_scan(sector, charge_filter, results_dir):
    model = EnergyModel()
    preset = preset_for_sector(sector)

    engine = ScanEngine(preset=preset, model=model)
    particles = get_particles()

    particles_filtered = {}
    for name, p in particles.items():
        if p.charge in charge_filter:
            particles_filtered[name] = p

    outputs = engine.run(particles_filtered)

    base_name = f"scan_{sector.value}_charge{''.join(charge_filter)}"
    png_path = results_dir / f"{base_name}.png"
    txt_path = results_dir / f"{base_name}.txt"
    csv_path = results_dir / f"{base_name}.csv"
    plot_scan(
        out_png=png_path,
        particles=particles_filtered,
        matched_points=outputs.matched_points,
        unmatched_segments=outputs.unmatched_segments,
        sector_name=sector.value,
    )

    write_results_txt(
        path=txt_path,
        sector=sector,
        particles=particles_filtered,
        bins_by_particle=outputs.bins_by_particle,
    )
    write_results_csv(
        path=csv_path,
        sector=sector,
        particles=particles_filtered,
        bins_by_particle=outputs.bins_by_particle,
    )

    png_grid_3d_path = results_dir / f"{base_name}_grid3d.png"
    plot_match_grid_3d_scatter(
        preset=preset, outputs=outputs, particles=particles_filtered, path=png_grid_3d_path
    )

    return png_path, png_grid_3d_path
