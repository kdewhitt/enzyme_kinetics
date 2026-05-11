from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure

_logger = logging.getLogger(__name__)
logging.getLogger("PIL").setLevel(logging.INFO)


def draw_kinetics_measured(
    ax: Axes,
    x: pd.Series,  # substrate_conc
    y: pd.Series,  # reaction_rate
    *,
    yerr: pd.Series,  # reaction_rate_err
    group_col: str,
    value_col: str = "intensity",
    config: Spec,
) -> None:

    _SCATTER_CONFIG = {
        "s": 5,  # Marker size
        "c": "#721482",  # Marker color
    }

    _SCATTER_ERROR_CONFIG = {
        "c": "#721482",  # Marker color
        "capsize": 2,
        "fmt": "r--o",
        "ecolor": "black",
        "linewidth": 1,
        "ms": 2,
    }

    ax.scatter(x, y, **_SCATTER_CONFIG)
    ax.errorbar(x, y, yerr=yerr, **_SCATTER_ERROR_CONFIG)

    if config.x_lims:
        ax.set_xlim(*config.x_lims)

    if config.y_lims:
        ax.set_ylim(*config.y_lims)

    if config.show_legend:
        ax.legend()


def draw_kinetics_model(
    ax: Axes,
    x: pd.Series,  # model_substrate_conc
    y: pd.Series,  # model_reaction_rate
    *,
    config: Spec,
) -> None:

    _LINE_CONFIG = {
        "color": "#138d49",
        "linewidth": 1,
    }

    ax.plot(x, y, **_LINE_CONFIG)

    if config.show_legend:
        ax.legend()


def make_kinetics_figure(
    substrate_conc: pd.Series,
    reaction_rate: pd.Series,
    reaction_rate_err: pd.Series,
    model_substrate_conc: pd.Series,
    model_reaction_rate: pd.Series,
    dest: str | Path | None = None,
    config: Spec | None = None,
) -> Figure:
    # Settings:
    # - width 3, height 2
    # - show_legend = True
    # - x_major, y_major = 10
    # - x_minor, y_minor = 2

    config = config or Spec()

    fig, ax = plt.subplots(figsize=config.figsize, layout="constrained")

    # Layer 3 - draw data
    draw_kinetics_measured(ax, substrate_conc, reaction_rate, yerr=reaction_rate_err, config=config)
    draw_kinetics_model(ax, model_substrate_conc, model_reaction_rate, config=config)

    # Layer 2b - apply style after data so tick locators see the data range
    apply_axes_style(ax, xlabel="Substrate concentration (mM)", ylabel="Velocity (mM/ min.)", config=config)

    # Layer 4 - save figure
    save_figure(fig, dest, description="kinetics scatter/line plot", config=config)

    return fig
