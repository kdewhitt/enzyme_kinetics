from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Final, Self

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from kgdlibs.pathtools import ExportPathBuilder, to_absolute_path
from plots.config import COUNT, MEAN, PEAK_ID, PROTEIN, REL_ACT, STD, filter_stats_params
from process.analyze import COLUMN_NAME_MAP, normalize_peak_id
from process.compounds import COMPOUND_PREFIXES
from process.io import read_csv_peaks
from process.normals import require_columns

from kinetics_core import (
    FitResult,
    KineticConstants,
    derive_constants,
    fit_michaelis_menten,
    fit_hill,
    fit_lineweaver_burk,
    michaelis_menten,
    hill_equation,
    prepare_velocity,
    fit_calibration,
    Calibration,
)

_logger = logging.getLogger(__name__)

SUBSTRATE_CONC: Final[str] = "substrate_conc"


# ---------------------------------------------------------------------------
# Config — data loading and export paths (unchanged from original)
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class KineticsConfig:
    prot_conc: float = 12.25
    rxn_time: float = 180.0
    peak_attr: str = "area"
    detect_prefix: bool = False
    peak_prefix: frozenset[str] = field(default_factory=frozenset)
    overwrite: bool = False
    merge: bool = True
    verbose: bool = False

    def __post_init__(self) -> None:
        prefixes = self.peak_prefix | frozenset(COMPOUND_PREFIXES)
        object.__setattr__(self, "peak_prefix", prefixes)

    def load(self, path: Path, **kwargs: Any) -> pd.DataFrame:
        df = read_csv_peaks(path, convert_to_float32=False, verbose=self.verbose, **kwargs)
        df.rename(columns=COLUMN_NAME_MAP, inplace=True)
        df = filter_stats_params(df, self.peak_attr)
        require_columns(df, (SUBSTRATE_CONC, PEAK_ID, MEAN, STD, COUNT, REL_ACT))
        prefixes = tuple(self.peak_prefix) if self.detect_prefix else None
        df = normalize_peak_id(df, prefix=prefixes)
        return df

    def build_export_path(self, dest: Path, tag: str) -> Path:
        base_name = f"{tag}_kinetics_{self.peak_attr}.png"
        build_path = ExportPathBuilder(
            base_name, overwrite=self.overwrite, merge=self.merge,
        ).build(dest)
        return build_path.path

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        valid_keys = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in valid_keys})


# ---------------------------------------------------------------------------
# Per-peak data extraction
# ---------------------------------------------------------------------------

def _extract_peak_data(
        df: pd.DataFrame,
        peak_id: str,
        rxn_time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    sub = df.loc[df[PEAK_ID] == peak_id].dropna(subset=[SUBSTRATE_CONC, MEAN])
    if sub.empty or sub[MEAN].max() == 0:
        return None
    return prepare_velocity(
        sub[SUBSTRATE_CONC].values,
        sub[MEAN].values,
        sub[STD].values,
        sub[COUNT].values,
        rxn_time,
    )


# ---------------------------------------------------------------------------
# Model selection strategy (replaces hardcoded 'OLV' branch)
# ---------------------------------------------------------------------------

def _select_and_fit(
        peak_id: str,
        s: np.ndarray,
        v: np.ndarray,
        v_sem: np.ndarray,
        *,
        special_peaks: dict[str, str] | None = None,
) -> FitResult:
    special = (special_peaks or {}).get(peak_id)
    if special == "hill":
        return fit_hill(s, v, sigma=v_sem)
    return fit_michaelis_menten(s, v, sigma=v_sem)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

class EnzymeKineticsAnalysis:
    def __init__(
            self,
            path: str | Path,
            target_dir: str | Path,
            *,
            config: KineticsConfig | None = None,
            special_peaks: dict[str, str] | None = None,
    ) -> None:
        self.path = to_absolute_path(path)
        self.dest = to_absolute_path(target_dir)
        self.config = config or KineticsConfig()
        self.special_peaks = special_peaks or {}

        self.df: pd.DataFrame = pd.DataFrame()
        self.results: dict[str, KineticConstants] = {}
        self.calibration: Calibration | None = None

    # -- pipeline stages ---------------------------------------------------

    def load(self, **kwargs: Any) -> Self:
        self.df = self.config.load(self.path, **kwargs)
        return self

    def fit(self) -> Self:
        self.results.clear()
        for peak_id in self.df[PEAK_ID].unique():
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            if data is None:
                continue
            s, v, v_sem = data
            try:
                mm_fit = _select_and_fit(
                    peak_id, s, v, v_sem, special_peaks=self.special_peaks,
                )
                lb_fit: FitResult | None = None
                try:
                    lb_fit = fit_lineweaver_burk(s, v)
                except Exception:
                    pass
                self.results[peak_id] = derive_constants(
                    peak_id, mm_fit, self.config.prot_conc, lb_fit=lb_fit,
                )
            except Exception as exc:
                _logger.warning("Fit failed for %s: %s", peak_id, exc)
        return self

    def apply_calibration(self, cal: Calibration) -> Self:
        self.calibration = cal
        recalculated: dict[str, KineticConstants] = {}
        for peak_id, kc in self.results.items():
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            if data is None:
                continue
            s, v, v_sem = data
            v_um = np.asarray(cal.area_to_conc(v * self.config.rxn_time)) / self.config.rxn_time
            v_sem_um = v_sem / cal.slope
            try:
                mm_fit = fit_michaelis_menten(s, v_um, sigma=v_sem_um)
                recalculated[peak_id] = derive_constants(
                    peak_id, mm_fit, self.config.prot_conc, lb_fit=kc.lb_fit,
                )
            except Exception as exc:
                _logger.warning("Calibrated re-fit failed for %s: %s", peak_id, exc)
        self.results = recalculated
        return self

    # -- output ------------------------------------------------------------

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([kc.to_dict() for kc in self.results.values()])

    def plot(self, *, show: bool = False, cols: int = 3) -> Self:
        peaks = list(self.results)
        if not peaks:
            _logger.warning("No fits to plot.")
            return self

        rows = int(np.ceil(len(peaks) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
        axes = np.atleast_1d(axes).flatten()

        for i, peak_id in enumerate(peaks):
            ax = axes[i]
            data = _extract_peak_data(self.df, peak_id, self.config.rxn_time)
            kc = self.results[peak_id]

            if data is None:
                ax.set_title(f"{peak_id} — no data")
                ax.axis("off")
                continue

            s, v, v_sem = data
            ax.errorbar(
                s, v, yerr=v_sem,
                fmt="o", color="steelblue", alpha=0.7, capsize=3, label="Data",
            )

            s_fit = np.linspace(0, np.max(s) * 1.1, 200)
            v_fit = kc.fit.predict(s_fit)
            label = (
                f"Vmax={kc.fit.vmax:.2f}±{kc.fit.vmax_se:.2f}\n"
                f"Km={kc.fit.km:.2f}±{kc.fit.km_se:.2f}\n"
                f"R²={kc.fit.r_squared:.4f}"
            )
            ax.plot(s_fit, v_fit, "--", color="red", label=label)
            ax.legend(fontsize=8, loc="lower right")
            ax.set_title(peak_id, fontweight="bold")
            ax.set_xlabel("[S] (µM)")
            ax.set_ylabel("v (signal / min)")

        for j in range(len(peaks), len(axes)):
            axes[j].axis("off")

        fig.tight_layout()
        plot_path = self.dest / f"{self.path.stem}_kinetics.png"
        fig.savefig(plot_path, dpi=300)
        if show:
            plt.show()
        plt.close(fig)
        _logger.info("Saved kinetics plot to %s", plot_path)
        return self

    def export_csv(self) -> Self:
        stats_df = self.to_dataframe()
        if stats_df.empty:
            _logger.warning("No kinetics data to export.")
            return self
        csv_path = self.config.build_export_path(
            self.dest, tag=self.path.stem,
        ).with_suffix(".csv")
        stats_df.to_csv(csv_path, index=False)
        _logger.info("Exported %s rows to %s", len(stats_df), csv_path)
        return self

    def run(self, **kwargs: Any) -> Self:
        return self.load(**kwargs).fit().plot().export_csv()
