from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import auto, IntFlag, StrEnum
from pathlib import Path
from typing import Any, Self

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from kgdlibs.datatools import map_columns
from kgdlibs.pathtools import ExportPathBuilder, to_absolute_path
from process.analyze_strip import COLUMN_NAME_MAP, normalize_peak_id
from process.infer import require_columns
from process.io import read_csv_peaks

_logger = logging.get_logger(__name__)


# TODO: 1ST PASS: Accept protein column
# TODO: 2ND PASS: Expect substrate and conc columns


class KineticsCols(StrEnum):
    SUBSTRATE = auto()
    CONC = auto()
    # PROTEIN = auto()
    PEAK_ID = auto()

    AREA_MEAN = auto()
    AREA_STD = auto()
    AREA_COUNT = auto()
    AREA_REL = auto()

    HEIGHT_MEAN = auto()
    HEIGHT_STD = auto()
    HEIGHT_COUNT = auto()
    HEIGHT_REL = auto()

    VMAX = auto()
    KM = auto()
    KM_UM = auto()
    KCAT = auto()
    KCAT_KM = auto()

    @classmethod
    def required_cols(cls) -> list[str]:
        """Returns the minimum set of columns required for analysis."""
        return [cls.SUBSTRATE, cls.CONC, cls.PEAK_ID]

    # protein?

    @classmethod
    def area_cols(cls) -> list[str]:
        return [cls.AREA_MEAN.value, cls.AREA_STD.value, cls.AREA_COUNT.value, cls.AREA_REL.value]

    @classmethod
    def height_cols(cls) -> list[str]:
        return [cls.HEIGHT_MEAN.value, cls.HEIGHT_STD.value, cls.HEIGHT_COUNT.value, cls.HEIGHT_REL.value]


# Column labels for readability
AREA = KineticsCols.AREA.value
HEIGHT = KineticsCols.HEIGHT.value
PEAK_ID = KineticsCols.PEAK_ID.value
# PROTEIN = KineticsCols.PROTEIN.value
# SAMPLE_ID = KineticsCols.SAMPLE_ID.value
# SAMPLE_NAME = KineticsCols.SAMPLE_NAME.value
SUBSTRATE = KineticsCols.SUBSTRATE.value


@dataclass(frozen=True, slots=True)
class KineticsConfig:
    protein_conc: float = 12.25

    # Import data
    case: str = "snake"
    scalar: float = 1.0
    convert_to_float32: bool = False
    include_cols: frozenset[str] = field(default_factory=lambda: frozenset((AREA, HEIGHT)))
    apply_filter: bool = False
    drop_na: bool = False
    verbose: bool = False

    # Normalize peak IDs
    detect_prefix: bool = False
    peak_prefix: frozenset[str] = field(default_factory=frozenset)

    overwrite: bool = False
    merge: bool = True

    @property
    def csv_args(self) -> dict[str, Any]:
        """A dictionary of keyword arguments forwarded to read_csv_peaks."""
        return {
            "case": self.case,
            "scalar": self.scalar,
            "convert_to_float32": self.convert_to_float32,
            "include_cols": self.include_cols,
            "apply_filter": self.apply_filter,
            "drop_na": self.drop_na,
            "verbose": self.verbose,
        }


def michaelis_menten(s, vmax, km):
    """Calculate the velocity term in the Michaelis–Menten kinetics equation."""
    return (vmax * s) / (km + s)
    # return (vmax * substrate_concentration) / (km + substrate_concentration)


def fit_kinetics(group, val_col: str = "mean", protein_conc: float = 1.0):
    s_data = group[SUBSTRATE].values
    v_data = group[val_col].values
    y_err = group["std"].values

    def loss(params):
        """Polynomial model loss calculation."""
        return np.sum((v_data - michaelis_menten(s_data, *params)) ** 2)

    # Initial guess [vmax, km]
    res = minimize(loss, [v_data.max(), np.median(s_data)], bounds=[(0, None), (1e-6, None)])
    # res = minimize(self._loss, [0, 1], product) where product is equivalent to group
    # self.s_model[product] = np.linspace(0, 4, 100)
    # self.v_model[product] = _calculate_velocity(self.s_model[product], res.x[0], res.x[1])

    vmax, km = res.x

    return pd.Series(
        {
            KineticsCols.VMAX.value: vmax,
            KineticsCols.KM.value: km,
            KineticsCols.KM_UM.value: km * 1000,
            KineticsCols.KCAT.value: vmax / (protein_conc / 10),
            KineticsCols.KCAT_KM.value: (vmax / (protein_conc / 10) / 60) / (km / 1000),
        },
    )


"""Calculate the kinetic parameters kcat and kcat/km.

            Variables:
                vmax = units of mM/min
                km = units of mM (per 3 hr)
                km_um = units of uM (per 1 hr)
                kcat = units of 1/min
                kcat/km = units of 1/s * 1/M
        """


class AreaHeightBitMask(IntFlag):
    NONE = 0
    AREA = 1
    HEIGHT = 2
    BOTH = AREA | HEIGHT


# self.s_real: dict[str, pd.Series] = {}  # {Product: Measured substrate concentrations}
# self.v_real: dict[str, pd.Series] = {}  # {Product: Calculated reaction velocity}
# self.y_err: dict[str, pd.Series] = {}  # {Product: Calculated error in y-values}
# self.s_model: dict[str, float | NDArray] = {}  # {Product: Modeled substrate concentrations}
# self.v_model: dict[str, float | NDArray] = {}  # {Product: Modeled reaction velocity}

class Kinetics:

    def __init__(
        self,
        path: str | Path,
        target_dir: str | Path,
        *,
        config: KineticsConfig | None = None,
    ) -> None:
        self.path = to_absolute_path(path)
        self.dest = to_absolute_path(target_dir)
        self.config = config or KineticsConfig()

        self.stats_df: pd.DataFrame = pd.DataFrame()

    def _load(self, **kwargs: Any):
        df = read_csv_peaks(self.path, **self.config.csv_args, **kwargs)
        df.rename(columns=COLUMN_NAME_MAP, inplace=True)
        # df['protein_clean'] = df['protein'].astype(str).str.replace(r"[a-zA-Z]+", "", regex=True).astype(float)
        # TODO: this was valid when protein column had values like "0.5mMMalCoA" --> Malcoa / 0.5 (separate cols)

        # Initial column validation
        # TODO: how to check for area_* and height_* columns?
        df = map_columns(df, KineticsCols.required_cols(), reindex=False)  # protein, peak_id, mean, std

        # Memory optimization and normalization
        prefixes = tuple(self.config.peak_prefix) if self.config.detect_prefix else None
        df = normalize_peak_id(df, prefix=prefixes)
        # df = normalize_protein_order(df, ref_protein=self.config.ref_protein)
        # df[SAMPLE_NAME] = df[SAMPLE_NAME].astype("category")
        # TODO: normalize substrate column

        return df

    def run(self) -> Self:
        df = self._load()

        bitmask = AreaHeightBitMask.NONE

        if not any(df.columns.str.contains("area") | df.columns.str.contains("height")):
            raise ValueError("DataFrame must contain area or height columns")

        if any(df.columns.str.contains("area")):
            bitmask += AreaHeightBitMask.AREA
            require_columns(df, set(KineticsCols.area_cols()))

        if any(df.columns.str.contains("height")):
            bitmask += AreaHeightBitMask.HEIGHT
            require_columns(df, set(KineticsCols.height_cols()))

        # Calculate stats using groupby-apply (highly optimized)
        self.stats_df = df.groupby(PEAK_ID).apply(
            fit_kinetics(protein_conc=self.config.protein_conc),
            include_groups=False,
        ).reset_index()

        return self

    def export(self) -> Self:
        if self.stats_df.empty:
            _logger.warning("PivotPeaks: No data to export.")
            return self

        build_path = ExportPathBuilder(
            f"{self.path.stem}_kinetics.csv",
            overwrite=self.config.overwrite,
            merge=self.config.merge,
        ).interpolate()

        self.stats_df.to_csv(build_path.path, index=False, index_label="index")

        _logger.info(
            "Saved enzyme kinetics data (shape=%s) to: %s",
            self.stats_df.shape, build_path.path,
        )

        return self
