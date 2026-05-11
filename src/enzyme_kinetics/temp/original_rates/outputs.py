import logging
from datetime import datetime
from pathlib import Path

from pandas import DataFrame

from config.constants import DATE_FORMAT
from old.kgdlibs_old import GetExportPath

_logger = logging.getLogger(__name__)


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

    def export_mm_statistics(self, data_frame: DataFrame) -> None:
        """Export the calculated Michaelis–Menten kinetic statistics to a CSV file."""
        file_name = f"{self.acquisition_date}_{self.path.stem}_kinetics_summary.csv"
        log_message = "Exported summarized Michaelis-Menten kinetics statistics to: \"{destination}\"."
        self._export_to_csv(data_frame, file_name, log_message)

    def _export_to_csv(self, data_frame: DataFrame, file_name: Path | str, log_message: str) -> Path:
        """General method to export data to CSV."""
        export_path_ = GetExportPath(self.export_path, file_name, self.overwrite, True)
        data_frame.to_csv(export_path_, index_label="index")
        _logger.info(log_message.format(destination=export_path_))
        return export_path_

    def suggest_plot_path(self, product: str) -> Path:
        """Determine an appropriate path to export a figure."""
        file_name = f"{self.acquisition_date}_{self.path.stem}_{product}_kinetics.png"
        return GetExportPath(self.export_path, file_name, self.overwrite, True)
