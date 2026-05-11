from enum import StrEnum


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
