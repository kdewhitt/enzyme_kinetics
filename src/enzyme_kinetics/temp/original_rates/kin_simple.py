from datetime import datetime
from pathlib import Path

from config import dict_to_class
from _messy_.kinetics.rates import CalculateEnzymeKinetics, ExportKineticsData, ImportKineticsData, PlotEnzymeKinetics


class KineticsHandler:
    """Handler for calculating and plotting enzyme kinetics data."""

    def __init__(
        self,
        path: Path | str,
        export_path: Path | str,
        date: datetime = None,
        overwrite: bool = False,
        show_figure: bool = False,
    ):
        """
        :param path: Path to a CSV file containing annotated HPLC peaks statistics to import and plot.
        :param export_path: Path to output processed data, figure(s), and/or log file(s).
        :param date: The data acquisition date usually in a date-like format, e.g. `20250104`.
        :param overwrite: Overwrite existing files.
        :param show_figure: Display figure in console.
        """
        self.path = Path(path).resolve()
        self.export_path = Path(export_path).resolve()
        self.run_date = date if date else datetime.today()
        self.overwrite = overwrite
        self.show_figure = show_figure

        # Class instance handlers
        self.h_data: ImportKineticsData | None = None
        self.h_calc: CalculateEnzymeKinetics | None = None
        self.h_export: ExportKineticsData | None = None

        # Perform operations
        self.h_export = ExportKineticsData(self.path, self.export_path, self.run_date, self.overwrite)

    def import_kinetics_data(self) -> None:
        """Import the annotated HPLC peaks statistics data."""
        self.h_data = ImportKineticsData(self.path).load()

    def calculate(self) -> None:
        """Calculate Michaelis-Menten enzyme kinetics."""
        self.h_calc = CalculateEnzymeKinetics(self.h_data.data, self.h_data.loc_names)
        self.h_calc.calculate_kinetics()
        self.h_calc.summarize()
        self.h_export.export_mm_statistics(self.h_calc.stats)

    def plot(self, params_dict: dict = None) -> None:
        """Plot Michaelis–Menten enzyme kinetics data as a line graph for each individual reaction product."""

        plotter = dict_to_class(PlotEnzymeKinetics, {**params_dict})

        for product in self.h_calc.unique:
            plotter.destination = self.h_export.suggest_plot_path(product)

            passed_arguments = {
                "substrate_conc": self.h_calc.s_real[product],
                "reaction_rate": self.h_calc.v_real[product],
                "reaction_rate_err": self.h_calc.y_err[product],
                "model_substrate_conc": self.h_calc.s_model[product],
                "model_reaction_rate": self.h_calc.v_model[product],
                "show_figure": self.show_figure,
            }

            plotter.kinetics_lineplot(**passed_arguments)
