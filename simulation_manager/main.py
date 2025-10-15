
import numpy as np
from physsm.rust_builder import RustExperimentBuilder, RustExperiment
from physsm.experiment_output import ExperimentOutput
from typing import override
from pathlib import Path

CARGO_TOML_PATH = Path.cwd().parent / "rust_simulation" / "Cargo.toml"
PROJECT_DIR     = Path.cwd().parent / "rust_simulation"


class IsingData(ExperimentOutput):
    def __init__(self, file_name):
        super().__init__(file_name)

        self.temperatures       = []
        self.energy_density     = []
        self.magnetisation      = []
        self.specific_heat      = []
        self.mag_susceptibility = []
        self.elapsed_time       = -1
        self.correlation_length = []

    @override
    def parse_output(self, line_number, line):
        if line_number == 0:
            slines = line.split(':')
            try:
                self.observables  = slines[0].split(', ')
                self.elapsed_time = float(slines[1])
            except Exception as _:
                self.observables = line.split(',')
                print("No elasped time found.")
        else:
            slines = line.split(", ")
            self.temperatures.append(float(slines[0]))
            self.energy_density.append(float(slines[1]))
            self.magnetisation.append(float(slines[2]))
            self.specific_heat.append(float(slines[3]))
            self.mag_susceptibility.append(float(slines[4]))
            self.correlation_length.append(float(slines[5]))
    
class RustIsingExperimentCreator:
    def __init__(self, folder: str, name: str):
        self.builder       = RustExperimentBuilder(PROJECT_DIR, folder, name)
        self.builder.set_scale_variable_names(["rows","cols"])
        self.builder.set_output_type(IsingData)

    def new_from_parameters(self, therm_steps: dict, measure_steps: dict, temperatures: np.ndarray, measure_struct_fact: bool = False) -> RustExperiment:
        
        self.builder.set_cargo_toml_path(CARGO_TOML_PATH)
        self.builder.add_static_parameter("temperatures", temperatures)
        self.builder.add_static_parameter("measure_struct_fact", measure_struct_fact)
        self.builder.add_scaling_parameter("therm_steps", therm_steps)
        self.builder.add_scaling_parameter("measure_steps", measure_steps)
        return self.builder.build()
    
    def load(self, lengths: list[int]) -> RustExperiment:
        self.builder.set_scale_variables(lengths)
        return self.builder.build(load_only=True)

def get_default_monte_carlo_parameters(lengths: list[int]) -> tuple[dict, dict]:
    """Parameters optimized for my machine (8Gb RAM, quadcore 2Ghz CPU (intel i7 8thGen))"""
    """Feel free to adapt for your machine"""
    
    thermalisation_steps = dict()
    measurement_steps    = dict()
    for L in lengths:
        if L <= 64:
            thermalisation_steps[L] = int(5e5)
            measurement_steps[L]    = int(5e5)
        elif L<= 128:
            thermalisation_steps[L] = int(1e5)
            measurement_steps[L]    = int(1e5)
        elif L<= 256:
            thermalisation_steps[L] = int(5e4) 
            measurement_steps[L]    = int(5e4)
        elif L<= 512:
            thermalisation_steps[L] = int(1e4) 
            measurement_steps[L]    = int(1e4)
        else:
            thermalisation_steps[L] = int(1e3) 
            measurement_steps[L]    = int(1e3)
    return (thermalisation_steps, measurement_steps)


def get_lengths(exp: RustExperiment) -> list[int]:
    return exp.get_scale_variables()

def get_results(exp: RustExperiment) -> dict[int, IsingData]:
    return exp.get_results()


def DEBUG_get_default_monte_carlo_parameters(lengths: list[int]) -> tuple[dict, dict]:
    """Debug params"""   
    thermalisation_steps = dict()
    measurement_steps    = dict()
    for L in lengths:
        thermalisation_steps[L] = int(1e3)
        measurement_steps[L]    = int(1e3)
    return (thermalisation_steps, measurement_steps)



if __name__ == "__main__":
    pass