import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pandas.core.arrays import ExtensionArray
from scipy.optimize import minimize

from config.constants import PROTEIN_CONCENTRATION
from old.kgdlibs_old.logtools import PPlog
from .constants import KineticsLabels, STATS_COLUMN_NAMES, StatsLabels


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
