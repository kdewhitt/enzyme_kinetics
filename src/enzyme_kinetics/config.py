from __future__ import annotations

from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Final, Self

import pandas as pd
from kgdlibs.pathtools import PathBuilder
from logurich import MultiLogAdapter

from enzyme_kinetics.compounds import canonicalize_peak_ids
from enzyme_kinetics.compounds.compounds import COMPOUND_PREFIXES

# _logger = logging.getLogger(__name__)
_logger = MultiLogAdapter(component="KineticsConfig", log_to_file=True)

# SUBSTRATE_CONC: Final[str] = "substrate_conc"
# COUNT = "count"
# MEAN = "mean"
# PEAK_ID = "peak_id"
# PROTEIN = "protein"
# REL_ACT = "rel_activity"
# STD = "std"


# ---------------------------------------------------------------------------
# Config — data loading and export paths
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class KineticsConfig:
    prot_conc: float = 12.25
    rxn_time: float = 180.0
    peak_attr: str = "area"
    detect_prefix: bool = False
    prefixes: frozenset[str] = field(default_factory=frozenset)
    overwrite: bool = False
    merge: bool = True
    verbose: bool = False

    def __post_init__(self) -> None:
        prefixes = self.prefixes | frozenset(COMPOUND_PREFIXES)
        object.__setattr__(self, "prefixes", prefixes)

    def load(self, path: Path, **kwargs: Any) -> pd.DataFrame:
        df = pd.read_csv(path, **kwargs)
        # df.rename(columns=COLUMN_NAME_MAP, inplace=True)
        # df = filter_stats_params(df, self.peak_attr)
        # require_columns(df, (SUBSTRATE_CONC, PEAK_ID, MEAN, STD, COUNT, REL_ACT))
        prefixes = tuple(self.prefixes) if self.detect_prefix else None
        df = canonicalize_peak_ids(df, prefixes=prefixes)
        return df

    def build_export_path(self, dest: Path, tag: str) -> Path:
        base_name = f"{tag}_kinetics_{self.peak_attr}.png"
        build_path = PathBuilder.from_path(
            dest / base_name,
            overwrite=self.config.overwrite,
        ).path
        return build_path

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Self:
        valid_keys = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in valid_keys})
