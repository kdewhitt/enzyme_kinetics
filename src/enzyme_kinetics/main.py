from pathlib import Path

from enzyme_kinetics.config import KineticsConfig


def example_kinetics_analysis():
    input_path_str = r"C:\Users\kris\Downloads\test_data_for_enzyme_kinetics_without_void.csv"
    input_path = Path(input_path_str)

    output_path = input_path.parent / "output"

    ENZYME_CONCENTRATION_UM = 12.25
    REACTION_TIME_HOURS = 3
    REACTION_TIME_SECONDS = REACTION_TIME_HOURS * 3600


    kinetics_config = KineticsConfig()
    loaded_df = kinetics_config.load(input_path)
    print(loaded_df.head())




if __name__ == "__main__":
    example_kinetics_analysis()
