import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import rcParams
from pandas import Series

from config import BaseFigure, modify_matplotlib_logging

_logger = logging.getLogger(__name__)

pil_logger = modify_matplotlib_logging()


class PlotEnzymeKinetics(BaseFigure):
    """Handler to plot Michaelis–Menten enzyme kinetics."""

    def __init__(self):
        super().__init__()

        self.x_label = "Substrate concentration (mM)"
        # self.x_label =  r"$[S]$ (mM)"
        self.y_label = "Velocity (mM/ min.)"
        # self.y_label = "Reaction rate (mM/min)"

        self.ticks.x_major = 10
        self.ticks.x_minor = 2
        self.ticks.y_major = 10
        self.ticks.y_minor = 2

    def kinetics_lineplot(
        self,
        substrate_conc: Series,
        reaction_rate: Series,
        reaction_rate_err: Series,
        model_substrate_conc: Series,
        model_reaction_rate: Series,
        show_figure: bool = False,
    ) -> None:
        """Plot Michaelis–Menten kinetics data as a line graph.

        :param substrate_conc: Substrate concentrations as x-values.
        :param reaction_rate: Reaction rates as y-values corresponding to substrate concentrations.
        :param reaction_rate_err: Errors for each reaction rate.
        :param model_substrate_conc: Substrate concentrations for the model.
        :param model_reaction_rate: Computed reaction rates for the model.
        :param show_figure: Display figure in console.
        """
        _apply_publication_matplotlib_global_style()

        # fig, ax = self.initialize_plot()
        fig, ax = plt.subplots(figsize=(3, 2), layout=self.figure.layout)

        # Plot raw data with error bar
        ax.scatter(substrate_conc, reaction_rate, c="#721482", s=5)
        # ax.scatter(substrate_conc, reaction_rate, marker='p', s=4)
        # ax.errorbar(substrate_conc, reaction_rate, yerr=reaction_rate_err,
        # capsize = 3, fmt = "r--o", ecolor = self.colors.error_bars)
        ax.errorbar(
            substrate_conc, reaction_rate, c="#721482", yerr=reaction_rate_err,
            capsize=2, fmt="r--o", ecolor=self.colors.error_bars, linewidth=1, ms=2,
        )

        # Plot model data
        ax.plot(model_substrate_conc, model_reaction_rate, color=self.colors.model_line, linewidth=1)

        # Modify figure, plot settings
        self.modify_axis_limits(ax)
        self.modify_axis_limits(ax)
        self.add_plot_labels(ax)

        ax.spines[['right', 'top']].set_visible(False)

        self.save_figure()
        self.save_figure(Path(self.destination).with_suffix(".eps"))
        if show_figure:
            plt.show()  # Show figure in console
        plt.clf()  # Clear matplotlib figure from memory


def _apply_publication_matplotlib_global_style() -> None:
    """Modify global matplotlib style settings."""
    rcParams["axes.grid"] = False
    rcParams["axes.xmargin"] = 0.00  # x-margin
    rcParams["font.family"] = "sans-serif"
    rcParams["font.sans-serif"] = "Arial"
    # rcParams["font.size"] = 5  # NEW!!
    rcParams["font.size"] = 6  # NEW!!
    # rcParams["axes.spines.left"] = False
    rcParams["savefig.dpi"] = 300  # Figure dots per inch or "figure"
