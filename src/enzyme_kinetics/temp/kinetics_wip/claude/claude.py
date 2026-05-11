"""Illustrative example of a layered matplotlib/seaborn plot architecture.

Demonstrates three-layer separation (style → layout → plot), a centralized
Spec dataclass, label mapping, and composable axis-level drawing
functions. Intended as a reference template — not runnable as-is.

Typical usage example:
    >>> config = Spec.publication()
    >>> fig, axes = make_two_panel_figure(config)
    >>> draw_intensity(axes[0], df, group_col="label", config=config)
    >>> draw_residuals(axes[1], df, config=config)
    >>> apply_axes_style(axes[0], xlabel="Time (s)", ylabel="Intensity", config=config)
    >>> apply_axes_style(axes[1], xlabel="Time (s)", ylabel="Residual", config=config)
    >>> fig.savefig("output.pdf")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# ---------------------------------------------------------------------------
# Label mapping — single source of truth for display strings
# ---------------------------------------------------------------------------

LABEL_MAP: Final[dict[str, str]] = {
    "HTAL": "H-Tal",
    "PDAL": "PD-Al",
    "OLA": "OL-A",
    "OLV": "OL-V",
    "TKS": "TKS",
    "M187A": "M187A",
    "M187A_ctrl": "M187A (ctrl)",
    "M187A_treated": "M187A (treated)",
}
"""Canonical display-label mapping for protein/sample identifiers."""


def _resolve_label(key: str) -> str:
    """Returns the display label for a protein or sample identifier."""
    return LABEL_MAP.get(key, key)


# ---------------------------------------------------------------------------
# Layer 1 — Style  (rcParams presets, palette helpers)
# ---------------------------------------------------------------------------

_RC_BASE: Final[dict] = {
    "font.family": "sans-serif",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}

_RC_PUBLICATION: Final[dict] = {
    **_RC_BASE,
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.linewidth": 0.6,
    "lines.linewidth": 0.9,
}

_RC_PRESENTATION: Final[dict] = {
    **_RC_BASE,
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
    "figure.dpi": 100,
}


def apply_rc_publication() -> None:
    """Applies publication-quality rcParams to the current matplotlib session.

    Configures tight spacing, small font sizes, and fine line weights suitable
    for journal figures. Must be called before figure creation to take effect.

    Examples:
        >>> apply_rc_publication()
        >>> fig, ax = plt.subplots()  # inherits publication settings
    """
    plt.rcParams.update(_RC_PUBLICATION)


def apply_rc_presentation() -> None:
    """Applies presentation-quality rcParams to the current matplotlib session.

    Configures larger fonts and heavier lines suitable for slides. Must be
    called before figure creation to take effect.

    Examples:
        >>> apply_rc_presentation()
        >>> fig, ax = plt.subplots()
    """
    plt.rcParams.update(_RC_PRESENTATION)


# ---------------------------------------------------------------------------
# Spec — declarative configuration passed through all layers
# ---------------------------------------------------------------------------

@dataclass(frozen=True, kw_only=True)
class Spec:
    """Immutable configuration bundle for a plot or set of related plots.

    Centralizes all cross-cutting options so drawing functions receive a
    single argument rather than a sprawling keyword list. Create instances
    via the named constructors (``publication()``, ``presentation()``) or
    by providing keyword arguments directly.

    Attributes:
        figsize: Figure dimensions in inches as (width, height).
        palette: Seaborn/matplotlib color palette name or list of colors.
        label_map: Mapping from raw data identifiers to display labels.
        x_major_interval: Spacing between major x-axis ticks.
        y_major_interval: Spacing between major y-axis ticks. If None,
            matplotlib chooses automatically.
        show_legend: If True, draws a legend on each axis that receives one.
        rc_preset: Name of the rcParams preset to apply on construction.
            One of ``"publication"``, ``"presentation"``, or ``"base"``.
    """

    figsize: tuple[float, float] = (6.5, 3.0)
    palette: str | list[str] = "colorblind"
    label_map: dict[str, str] = field(default_factory=lambda: dict(LABEL_MAP))
    x_major_interval: float | None = None
    y_major_interval: float | None = None
    show_legend: bool = True
    rc_preset: str = "publication"

    def __post_init__(self) -> None:
        """Applies the configured rcParams preset after initialization."""
        match self.rc_preset:
            case "publication":
                apply_rc_publication()
            case "presentation":
                apply_rc_presentation()
            case "base":
                plt.rcParams.update(_RC_BASE)
            case _:
                raise ValueError(
                    f"Unknown rc_preset {self.rc_preset!r}. "
                    "Expected 'publication', 'presentation', or 'base'.",
                )

    def resolve_label(self, key: str) -> str:
        """Returns the display label for a data identifier."""
        return self.label_map.get(key, key)

    @classmethod
    def publication(cls, **overrides) -> Spec:
        """Constructs a config preset tuned for journal figures.

        Args:
            **overrides: Any Spec field values to override.

        Returns:
            A Spec instance with publication defaults applied.

        Examples:
            >>> cfg = Spec.publication(figsize=(3.5, 2.5))
        """
        return cls(rc_preset="publication", figsize=(6.5, 3.0), **overrides)

    @classmethod
    def presentation(cls, **overrides) -> Spec:
        """Constructs a config preset tuned for slide figures.

        Args:
            **overrides: Any Spec field values to override.

        Returns:
            A Spec instance with presentation defaults applied.

        Examples:
            >>> cfg = Spec.presentation(show_legend=False)
        """
        return cls(rc_preset="presentation", figsize=(10.0, 5.0), **overrides)


# ---------------------------------------------------------------------------
# Layer 2 — Layout  (figure/axes factories, no data)
# ---------------------------------------------------------------------------

def make_figure(
    config: Spec,
    *,
    nrows: int = 1,
    ncols: int = 1,
    sharex: bool = False,
    sharey: bool = False,
) -> tuple[Figure, list[Axes]]:
    """Creates a figure with a grid of axes using the provided config.

    Returns axes as a flat list regardless of grid shape, which simplifies
    iteration over panels without needing to handle 1-D vs 2-D numpy arrays.

    Args:
        config: Plot configuration supplying figsize and rcParams preset.
        nrows: Number of subplot rows.
        ncols: Number of subplot columns.
        sharex: If True, all axes share the same x-axis limits and ticks.
        sharey: If True, all axes share the same y-axis limits and ticks.

    Returns:
        A tuple of (Figure, flat list of Axes) with nrows * ncols elements.

    Examples:
        >>> cfg = Spec.publication()
        >>> fig, axes = make_figure(cfg, nrows=1, ncols=2, sharex=True)
        >>> len(axes)
        2
    """
    fig, ax_array = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=config.figsize,
        sharex=sharex,
        sharey=sharey,
        constrained_layout=True,
    )
    axes = [ax_array] if nrows == 1 and ncols == 1 else list(ax_array.flat)
    return fig, axes


# ---------------------------------------------------------------------------
# Layer 2b — Axes style  (aesthetics applied after layout, before or after data)
# ---------------------------------------------------------------------------

def apply_axes_style(
    ax: Axes,
    *,
    xlabel: str = "",
    ylabel: str = "",
    title: str = "",
    config: Spec,
) -> None:
    """Applies consistent axis-level formatting to a single Axes object.

    Configures tick intervals, spine visibility, and axis labels using values
    from the provided config. Separates aesthetic concerns from data drawing
    so plot functions stay focused on rendering data.

    Args:
        ax: The Axes instance to format.
        xlabel: Label for the x-axis.
        ylabel: Label for the y-axis.
        title: Axes title text. Use an empty string to suppress.
        config: Plot configuration supplying tick interval and style options.

    Examples:
        >>> apply_axes_style(ax, xlabel="Time (s)", ylabel="AU", config=cfg)
    """
    if config.x_major_interval is not None:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(config.x_major_interval))
    if config.y_major_interval is not None:
        ax.yaxis.set_major_locator(ticker.MultipleLocator(config.y_major_interval))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    # Remove top/right spines for a cleaner publication look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


# ---------------------------------------------------------------------------
# Layer 3 — Plot  (axis-level drawing functions, data in + axes in)
# ---------------------------------------------------------------------------

def draw_intensity(
    ax: Axes,
    data,  # pd.DataFrame
    *,
    group_col: str,
    value_col: str = "intensity",
    config: Spec,
) -> None:
    """Draws a line plot of intensity traces grouped by sample label.

    Each group is drawn as a separate line, with labels resolved through
    config.label_map. Legend is drawn when config.show_legend is True.
    Never calls plt.show() or modifies rcParams — all state is local to ax.

    Args:
        ax: Target Axes to draw onto.
        data: DataFrame containing at minimum group_col, value_col, and a
            column named "time" for the x-axis.
        group_col: Column name identifying the sample or protein group.
        value_col: Column name for the y-axis values. Defaults to "intensity".
        config: Plot configuration supplying palette and label options.

    Examples:
        >>> draw_intensity(ax, df, group_col="label", config=cfg)
    """
    groups = data[group_col].unique()
    colors = sns.color_palette(config.palette, n_colors=len(groups))

    for group, color in zip(groups, colors):
        subset = data[data[group_col] == group]
        label = config.resolve_label(group)
        ax.plot(subset["time"], subset[value_col], label=label, color=color)

    if config.show_legend:
        ax.legend()


def draw_bar_comparison(
    ax: Axes,
    data,  # pd.DataFrame
    *,
    group_col: str,
    value_col: str,
    error_col: str | None = None,
    config: Spec,
) -> None:
    """Draws a bar chart comparing values across groups with optional error bar.

    Group labels are resolved through config.label_map before rendering.
    The function draws onto the provided axes without modifying figure-level
    state or calling plt.show().

    Args:
        ax: Target Axes to draw onto.
        data: DataFrame with one row per group, or a summary DataFrame with
            pre-computed means and errors.
        group_col: Column identifying each group.
        value_col: Column holding the bar heights.
        error_col: Column holding symmetric error bar half-widths. If None,
            error bar are omitted.
        config: Plot configuration supplying palette and label options.

    Examples:
        >>> draw_bar_comparison(ax, summary_df, group_col="sample",
        ...                     value_col="mean_intensity", config=cfg)
    """
    labels = [config.resolve_label(g) for g in data[group_col]]
    colors = sns.color_palette(config.palette, n_colors=len(labels))
    yerr = data[error_col].to_numpy() if error_col else None

    ax.bar(
        labels,
        data[value_col],
        yerr=yerr,
        color=colors,
        capsize=3,
        error_kw={"linewidth": 0.8},
    )


# ---------------------------------------------------------------------------
# Composition example — assembles layers into a complete figure
# ---------------------------------------------------------------------------

def make_intensity_bar_figure(
    intensity_df,  # pd.DataFrame
    summary_df,  # pd.DataFrame
    *,
    group_col: str = "label",
    config: Spec | None = None,
) -> Figure:
    """Produces a two-panel figure with intensity traces and a bar comparison.

    Demonstrates full layer composition: layout factory → drawing functions →
    axes style. The caller receives the Figure and can save or further annotate
    it. rcParams are applied via config before figure creation.

    Args:
        intensity_df: Long-form DataFrame with time-series intensity data.
        summary_df: Summary DataFrame with pre-computed group means and errors.
        group_col: Column name identifying sample/protein groups in both
            DataFrames. Defaults to "label".
        config: Plot configuration to use. If None, uses
            ``Spec.publication()`` with default settings.

    Returns:
        The completed Figure with two panels ready for saving.

    Examples:
        >>> fig = make_intensity_bar_figure(
        ...     intensity_df, summary_df, config=Spec.publication()
        ... )
        >>> fig.savefig("comparison.pdf")
    """
    if config is None:
        config = Spec.publication()

    fig, (ax_trace, ax_bar) = plt.subplots(
        1, 2, figsize=config.figsize, constrained_layout=True,
    )

    # Layer 3 — draw data
    draw_intensity(ax_trace, intensity_df, group_col=group_col, config=config)
    draw_bar_comparison(
        ax_bar,
        summary_df,
        group_col=group_col,
        value_col="mean",
        error_col="sem",
        config=config,
    )

    # Layer 2b — apply style after data so tick locators see the data range
    apply_axes_style(
        ax_trace,
        xlabel="Time (s)",
        ylabel="Intensity (AU)",
        title="Traces",
        config=config,
    )
    apply_axes_style(
        ax_bar,
        ylabel="Mean Intensity (AU)",
        title="Comparison",
        config=config,
    )

    return fig
