import logging
from datetime import datetime
from enum import StrEnum
from pathlib import Path

import numpy as np
import pandas as pd
from config.constants import CASE_STYLE, DATE_FORMAT, DROP_NA, PROTEIN_CONCENTRATION
from numpy.typing import NDArray
from pandas.core.arrays import ExtensionArray
from scipy.optimize import minimize

from with_docs.kgdlibs import GetExportPath
from with_docs.kgdlibs.datatools import describe_data, ImportCSVFile, MapColumnNamesByCase
from with_docs.kgdlibs.logtools import PPlog

_logger = logging.getLogger(__name__)


class KineticsLabels(StrEnum):
    # INDEX = "index"
    PROTEIN = "protein"
    SUBSTRATE = "substrate"
    PEAK_ID = "peak_id"
    MEAN = "mean"
    STD = "std"
    COUNT = "count"
    PERCENT_ACTIVITY = "percent_activity"


ALL_COLUMN_NAMES = [
    # KineticsLabels.INDEX,
    KineticsLabels.PROTEIN,
    KineticsLabels.SUBSTRATE,
    KineticsLabels.PEAK_ID,
    KineticsLabels.MEAN,
    KineticsLabels.STD,
    KineticsLabels.COUNT,
    KineticsLabels.PERCENT_ACTIVITY,
]


class StatsLabels(StrEnum):
    PRODUCT = "product"
    VMAX = "vmax"
    KM = "km"
    KM_UM = "km_um"
    KCAT = "kcat"
    KCAT_KM = "kcat_km"


STATS_COLUMN_NAMES = [
    StatsLabels.PRODUCT,
    StatsLabels.VMAX,
    StatsLabels.KM,
]

_ALPHA_PATTERN = r"[a-zA-Z]+"


class ImportKineticsData:
    """Import annotated HPLC peaks statistics."""

    def __init__(self, path: Path | str):
        """
        :param path: Path to a CSV file containing annotated HPLC peaks statistics to import and plot.
        """
        self.path = Path(path).resolve()

        self.data: pd.DataFrame = pd.DataFrame()  # Imported data
        self.loc_names: dict = {}  # Mapped column names // ALL_COLUMN_NAMES

    def load(self):
        """Import annotated HPLC peaks statistics, standardize column names, and preprocess the data."""
        self.data = self._load()
        self.loc_names = self._map_column_names(ALL_COLUMN_NAMES)
        self._clean_protein_col()
        self._describe_peaks_data()
        return self

    def _load(self) -> pd.DataFrame:
        """Import HPLC peaks data and standardize column names."""
        return ImportCSVFile(CASE_STYLE, DROP_NA).load(self.path)

    def _map_column_names(self, col_names: list[str]) -> dict:
        """Map reference column names to case-converted column names."""
        return MapColumnNamesByCase(col_names, CASE_STYLE).perform_mapping(self.data)  # Mapped names

    def _clean_protein_col(self):
        """
        Clean the column labeled `protein` by removing alphabetic
        characters and converting the remaining values to floats.
        """
        self.data.loc[:, self.loc_names[KineticsLabels.PROTEIN]] = (
            _strip_alpha_chars(self.data, KineticsLabels.PROTEIN)
        )

    def _describe_peaks_data(self) -> pd.DataFrame:
        """Describe the imported data."""
        describe_data(self.data)
        return self.data


def _strip_alpha_chars(data_frame: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """Remove alphabetic characters and converting the remaining values to floats."""
    # data_frame.loc[:, "protein"] = data_frame.loc[:, col_name].str.replace(_ALPHA_PATTERN, "", regex=True).astype(
    # float)
    return data_frame.loc[:, col_name].str.replace(_ALPHA_PATTERN, "", regex=True).astype(float)


class CalculateEnzymeKinetics:
    """Calculate Michaelis–Menten enzyme kinetics."""

    def __init__(self, imported_data: pd.DataFrame, mapping_dict: dict):
        """
        :param imported_data: The imported peaks data to use for calculating peak statistics.
        :param mapping_dict: Map data frame columns to case-converted names.
        """
        self.data: pd.DataFrame = imported_data
        self.loc_names: dict = mapping_dict  # Mapped column names // ALL_COLUMN_NAMES

        self.unique: ExtensionArray | NDArray = np.array([])  # Unique products

        self.s_real: dict[str, pd.Series] = {}  # {Product: Measured substrate concentrations}
        self.v_real: dict[str, pd.Series] = {}  # {Product: Calculated reaction velocity}
        self.y_err: dict[str, pd.Series] = {}  # {Product: Calculated error in y-values}
        self.s_model: dict[str, float | NDArray] = {}  # {Product: Modeled substrate concentrations}
        self.v_model: dict[str, float | NDArray] = {}  # {Product: Modeled reaction velocity}

        self.stats: pd.DataFrame = pd.DataFrame(
            columns=STATS_COLUMN_NAMES,
        )  # Michaelis–Menten statistics, e.g. vmax, km

        # Perform operations
        self._get_unique_products()

    # ///////////////////
    # ///////////////////
    # ///////////////////
    # ///////////////////
    # 20260510
    # 20260510
    # 20260510
    # 20260510
    # ///////////////////
    # ///////////////////
    # ///////////////////
    # ///////////////////


    def _get_unique_products(self) -> None:
        """Return a unique list of products, e.g. HTAL, PDAL, olivetolic acid, olivetol."""
        self.unique = self.data[self.loc_names[KineticsLabels.PEAK_ID]].unique()

    def calculate_kinetics(self) -> None:
        """Calculate Michaelis–Menten kinetics for each unique product."""
        for product in self.unique:
            product_data = self.data[self.data[self.loc_names[KineticsLabels.PEAK_ID]] == product]
            self.s_real[product] = product_data.loc[:, self.loc_names[KineticsLabels.PROTEIN]]
            self.v_real[product] = product_data.loc[:, self.loc_names[KineticsLabels.MEAN]]
            # self.v_real[product] = product_data.loc[:, self.loc_names[KineticsLabels.PERCENT_ACTIVITY]]
            self.y_err[product] = product_data.loc[:, self.loc_names[KineticsLabels.STD]]

            res = minimize(self._loss, [0, 1], product)
            self.s_model[product] = np.linspace(0, 4, 100)
            self.v_model[product] = _calculate_velocity(self.s_model[product], res.x[0], res.x[1])

            # # self.stats = self.stats.append(
            # self.stats = pd.concat([self.stats,
            #     {
            #         StatsLabels.PRODUCT: product,
            #         StatsLabels.VMAX: res.x[0],
            #         StatsLabels.KM: res.x[1]
            #     }],axis=1,
            #     ignore_index=True)

            dat1 = pd.DataFrame(
                {
                    StatsLabels.PRODUCT: product,
                    StatsLabels.VMAX: res.x[0],
                    StatsLabels.KM: res.x[1],
                }, index=[0],
            )

            print(dat1)
            self.stats = pd.concat([self.stats, dat1])
            print(self.stats)

            PPlog.debug(f"Calculated error bars for \"{product}\" are:", self.y_err[product])
            PPlog.debug(f"Fitted model's x-values for \"{product}\" are:", res.x)

    def summarize(self) -> None:
        """Calculate the kinetic parameters kcat and kcat/km.

            Variables:
                vmax = units of mM/min
                km = units of mM (per 3 hr)
                km_um = units of uM (per 1 hr)
                kcat = units of 1/min
                kcat/km = units of 1/s * 1/M
        """
        self.stats[StatsLabels.KM_UM] = self._calculate_km_um()
        self.stats[StatsLabels.KCAT] = self._calculate_kcat()
        self.stats[StatsLabels.KCAT_KM] = self._calculate_kcat_km()

    def _loss(self, params, product: str) -> float:
        """Polynomial model loss calculation."""
        v_max, k_m = params
        v_pred = _calculate_velocity(self.s_real[product], v_max, k_m)
        return np.sum((self.v_real[product] - v_pred) ** 2)

    def _calculate_km_um(self) -> float | pd.Series | NDArray | ExtensionArray | None:
        """Convert Km units from mM to uM."""
        return self.stats.loc[:, "km"] * 1000

    def _calculate_kcat(self) -> float | pd.Series | NDArray | ExtensionArray | None:
        """Calculate kcat from vmax."""
        return self.stats.loc[:, "vmax"] / (PROTEIN_CONCENTRATION / 10)

    def _calculate_kcat_km(self) -> float | pd.Series | NDArray | ExtensionArray | None:
        """Calculate kcat/km."""
        return (self.stats.loc[:, "kcat"] / 60) / (self.stats["km"] / 1000)


def _calculate_velocity(
    substrate_concentration: int | float | pd.Series,
    vmax: int | float, km: int | float,
) -> float:
    """Calculate the velocity term in the Michaelis–Menten kinetics equation."""
    return (vmax * substrate_concentration) / (km + substrate_concentration)


class ExportKineticsData:

    def __init__(
        self,
        path: Path | str,
        export_path: Path | str,
        run_date: datetime,
        overwrite: bool = False,
    ):
        """
        :param path: Path to a CSV file containing annotated HPLC peaks statistics to export.
        :param export_path: Path to output processed data, figure(s), and/or log file(s).
        :param run_date: The data acquisition date usually in a date-like format, e.g. `20250104`.
        :param overwrite: Overwrite existing files.
        """
        self.path = Path(path).resolve()
        self.export_path = Path(export_path).resolve()
        self.run_date = run_date
        self.overwrite = overwrite

    @property
    def acquisition_date(self) -> str:
        """Return the data acquisition date (i.e. `run_date`) formatted as `YYYYMMDD`."""
        return self.run_date.strftime(DATE_FORMAT)

    def export_mm_statistics(self, data_frame: pd.DataFrame) -> None:
        """Export the calculated Michaelis–Menten kinetic statistics to a CSV file."""
        file_name = f"{self.acquisition_date}_{self.path.stem}_kinetics_summary.csv"
        log_message = "Exported summarized Michaelis-Menten kinetics statistics to: \"{destination}\"."
        self._export_to_csv(data_frame, file_name, log_message)

    def _export_to_csv(self, data_frame: pd.DataFrame, file_name: Path | str, log_message: str) -> Path:
        """General method to export data to CSV."""
        export_path_ = GetExportPath(self.export_path, file_name, self.overwrite, True)
        data_frame.to_csv(export_path_, index_label="index")
        _logger.info(log_message.format(destination=export_path_))
        return export_path_

    def suggest_plot_path(self, product: str) -> Path:
        """Determine an appropriate path to export a figure."""
        file_name = f"{self.acquisition_date}_{self.path.stem}_{product}_kinetics.png"
        return GetExportPath(self.export_path, file_name, self.overwrite, True)
