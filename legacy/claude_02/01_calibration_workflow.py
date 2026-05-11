#!/usr/bin/env python3
"""Example: CoA Calibration Workflow

Demonstrates how to:
1. Load calibration standards
2. Fit calibration curve
3. Validate fit quality
4. Save calibration for reuse
"""

from pathlib import Path
import numpy as np

# Add parent directory to path for imports
import sys

sys.path.insert(0, str(Path(__file__).parent))

from enzyme_kinetics import CalibrationCurve


def example_calibration_workflow():
    """Demonstrates complete calibration workflow."""
    print("=" * 80)
    print("CALIBRATION WORKFLOW EXAMPLE")
    print("=" * 80 + "\n")

    # Step 1: Load or create calibration data
    # In practice, this would come from your HPLC measurements
    print("Step 1: Loading calibration standards...")

    # Example calibration data (replace with your actual data)
    calibration_data = {
        "coa_concentration_um": np.array([0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]),
        "peak_area_mean": np.array([0, 245.3, 612.4, 1223.1, 2451.8, 4903.5, 12258.7]),
        "peak_area_std": np.array([0, 12.4, 28.1, 45.2, 89.3, 156.2, 412.1]),
    }

    conc = calibration_data["coa_concentration_um"]
    area = calibration_data["peak_area_mean"]
    area_std = calibration_data["peak_area_std"]

    print(f"  Loaded {len(conc)} calibration points")
    print(f"  Concentration range: {conc.min():.2f} to {conc.max():.2f} µM\n")

    # Step 2: Fit calibration curve
    print("Step 2: Fitting calibration curve...")

    calib = CalibrationCurve()
    params = calib.fit(conc, area, peak_areas_std=area_std)

    print(f"  {params}\n")

    # Step 3: Validate fit
    print("Step 3: Validating fit quality...")

    if calib.is_valid():
        print(f"  ✓ Fit is valid (R² = {params.r2:.6f})\n")
    else:
        print(f"  ⚠ Fit may not be reliable (R² = {params.r2:.6f})\n")

    # Step 4: Test conversion functions
    print("Step 4: Testing area-to-concentration conversion...")

    test_areas = np.array([245.3, 612.4, 2451.8])
    test_conc = calib.area_to_concentration(test_areas)

    print(f"  Area 245.3 → Concentration {test_conc[0]:.4f} µM")
    print(f"  Area 612.4 → Concentration {test_conc[1]:.4f} µM")
    print(f"  Area 2451.8 → Concentration {test_conc[2]:.4f} µM\n")

    # Step 5: Error propagation example
    print("Step 5: Testing error propagation...")

    area_measurement = 612.4
    area_std_measurement = 28.1
    conc, conc_std = calib.area_to_concentration_with_error(
        area_measurement, area_std_measurement
    )

    print(f"  Area: {area_measurement} ± {area_std_measurement}")
    print(f"  Concentration: {conc:.4f} ± {conc_std:.4f} µM\n")

    # Step 6: Save calibration
    print("Step 6: Saving calibration parameters...")

    output_dir = Path("/mnt/user-data/outputs")
    calib_file = output_dir / "calibration_example.txt"
    calib.save(calib_file)
    print(f"  ✓ Saved to {calib_file}\n")

    # Step 7: Generate plots
    print("Step 7: Generating calibration plots...")

    plot_path = output_dir / "calibration_example.png"
    calib.plot(output_path=plot_path)
    print()

    print("=" * 80)
    print("✓ CALIBRATION WORKFLOW COMPLETE")
    print("=" * 80 + "\n")

    return calib


if __name__ == "__main__":
    calib = example_calibration_workflow()
