from __future__ import annotations
import numpy as np
import os
import time
from .io_handler import perform_rust_computation, write_parameter_file
from .rust_montecarlo_data_handler import RustMonteCarloData

Parameter = dict[int, int]

def array_to_str(array: np.ndarray, rounding: int) -> list[str]:
    temps_str = []
    for val in array:
        temps_str.append(str(round(val,rounding)))
    return ", ".join(temps_str)

class RustExperiment:
    
    def __init__(self, name: str, folder: str, lengths):
        self.name: str                  = name
        self.folder: str                = folder
        self.lengths: np.ndarray        = lengths
        self.measure_struct_factor:bool = False
        self.therm_steps: np.ndarray    = None
        self.measure_steps: np.ndarray  = None
        self.temperatures: np.ndarray   = None
        self.param_files: dict          = dict()
        self.out_files: dict            = dict()
        self.init_directory()
        self.set_files()
    
    def write_parameter_file(self, L: int, rounding = 3) -> None:
        if self.temperatures is None:
            print("No temperatures set!")
            return
        if self.therm_steps is None or self.measure_steps is None:
            print("Monte Carlo Parameters not set!")
            return
        if self.param_files is None or self.out_files is None:
            print("Parameter & output files not set!")
        try:
            with open(self.param_files[L], "w") as f:
                f.write(f"rows: {L}\n")
                f.write(f"cols: {L}\n")
                f.write(f"therm_steps: {self.therm_steps[L]}\n")
                f.write(f"measure_steps: {self.measure_steps[L]}\n")
                f.write(f"temperatures: {array_to_str(self.temperatures, rounding)}\n")
                f.write(f"measure_struct_fact: {self.measure_struct_factor}\n")
                f.write(f"outputfile: {self.out_files[L]}\n")
            print(f"-- {self.param_files[L]}: ", end="")
        except Exception as e:
            print(f"Error occured: {e}, {e.args}")
                
    def write_parameter_files(self, temp_rounding=3):
        print(f">> writing parameter files:")
        for L in self.lengths:
            try:
                self.write_parameter_file(L, temp_rounding)
                print("done")
            except Exception as e:
                print(f"Error occured: {e}, {e.args}")
                
    def set_files(self):
        print(">> Setting file paths:")
        if self.lengths is None:
            print(">> Set Files: Error: No lengths found!")
            return

        for L in self.lengths:
            self.param_files[L] = f"{self.folder}/{self.name}/parameter_{L}x{L}.txt"
            print(f"-- setting paramfile \"{self.param_files[L]}\"")
        for L in self.lengths:
            self.out_files[L]   = f"{self.folder}/{self.name}/out_{L}x{L}.txt"
            print(f"-- setting outfile \"{self.out_files[L]}\"")
    
    def should_measure_struct_fact(self, value: bool):
        self.measure_struct_factor = value
        print(f">> measuring structur factor: {value}")
    
    def has_param_file(self, L: int) -> bool:
        try:
            with open(f"{self.folder}/{self.name}/parameter_{L}x{L}.txt", "r") as _:
                pass
            return True
        except FileNotFoundError as _:
            return False
    
    def are_parameter_files_available(self) -> bool:
        has_all_files = True
        missing_files = []
        print(">> Looking for parameter files:")
        if self.lengths is None:
            print(f"* No lengths set")
            return False
        
        for L in self.lengths:
            file = f"\"{self.folder}/{self.name}/parameter_{L}x{L}.txt\""
            if self.has_param_file(L):
                print(f"-- Found parameter file: {file}")
            else:
                missing_files.append(file)
                has_all_files = False
                
        if has_all_files:
            print("* Parameter files available")
        else:
            print("* Not all parameter files created yet:")    
            for files in missing_files:
                print(f"-- not written yet: {files}")
        return has_all_files

    
    @staticmethod
    def new_from_parameters(name: str, folder: str, therm_steps: Parameter, measure_steps: Parameter, temperatures: np.ndarray):
        lengths = RustExperiment.get_lengths_from_parameters(therm_steps, measure_steps)
        
        new_exp = RustExperiment(name, folder, lengths)
        new_exp.set_parameters(therm_steps, measure_steps)
        new_exp.set_temperatures(temperatures)
        return new_exp
    
    @staticmethod
    def get_lengths_from_parameters(therm_steps: Parameter, measure_steps: Parameter):
        lengths = []
        for (L, _) in zip(therm_steps.keys(), measure_steps.keys()):
            lengths.append(L)
        return lengths
        
    def set_parameters(self, therm_steps: Parameter, measure_steps: Parameter):
        param_lengths = RustExperiment.get_lengths_from_parameters(therm_steps, measure_steps)
        if param_lengths != self.lengths:
            print(">> Monte carlo parameter lengths do not match original lengths!")
        else:
            self.measure_steps = measure_steps
            self.therm_steps   = therm_steps
        print(">> Monte Carlo parameters set!")
        
    # def set_parameters_default(self):
    #     for L in self.lengths:
    #         if L <= 64:
    #             self.therm_steps[L] = np.uint(5e5)
    #             self.measure_steps[L]    = np.uint(5e5)
    #         elif L <= 128:
    #             self.therm_steps[L] = np.uint(1e5)
    #             self.measure_steps[L]    = np.uint(1e5)
    #         else:
    #             self.therm_steps[L] = np.uint(1e4)
    #             self.measure_steps[L]    = np.uint(1e4)
    
    def are_results_available(self):
        all_available = True
        some_availabe = False
        missing_files = []
        print(">> Looking for output files:")
        for L in self.lengths:
            file = f"\"{self.folder}/{self.name}/out_{L}x{L}.txt\""
            if self.has_output(L):
                print(f"-- Found output: {file}")
                some_availabe = True
            else:
                missing_files.append(file)
                all_available = False
                
        if not all_available:
            print(f"* missing output:")
            for files in missing_files:
                print(f"-- {files}")
        if not some_availabe:
            print("* No output files created yet.")
            
        return some_availabe
        
    def init_directory(self):
        try:
            os.makedirs(f"{self.folder}/{self.name}")
            print(f">> Directory \"{self.folder}/{self.name}\" created")
        except FileExistsError:
            print(f">> Directory \"{self.folder}/{self.name}\" Found.")
        except OSError as e:
            print(f"Error creating directory: {e}")

            
    def set_temperatures(self, temps: np.ndarray):
        if isinstance(temps, np.ndarray):
            print(">> Temperatures set!")
            self.temperatures = temps
        else:
            print(">> Temperatures not numpy array!")
    # def write_parameter_file(self, L):
    #     if self.temperatures is None:
    #             print("No temperatures set!")
    #             return
    #         if self.therm_steps is None or self.measure_steps is None:
    #             print("Monte Carlo Parameters not set!")
    #             return
    #         if self.param_files is None or self.out_files is None:
    #             print("Parameter & output files not set!")
    #     try:
    #         print(f">> writing: {self.param_files[L]}")
    #         write_parameter_file(self.param_files[L], self.out_files[L], L, L, self.therm_steps[L], self.measure_steps[L], self.temperatures)
    #     except Exception as e:
    #         print(f"Error occured: {e}, {e.args}")
                

                

    def perform_rust_calculation(self, L: int) -> int:
        'Returns time in secs'        
        
        if L not in self.lengths:
            print(f"Length {L} not in experiment")
        return perform_rust_computation(self.param_files[L])
    
    def has_output(self, L):
        try:
            with open(f"{self.folder}/{self.name}/out_{L}x{L}.txt", "r") as _:
                pass
            # print(f"Found output: \"{self.folder}/{self.name}/out_{L}x{L}.txt")
            return True
        except FileNotFoundError as e:
            return False
    

    
    def get_results(self) -> dict[int, RustMonteCarloData]: #tuple[np.ndarray, dict[int, RustMonteCarloData]]:
        print("")
        results: dict[int, RustMonteCarloData] = dict()
        temp_set = False
        for L in self.lengths:
            if self.has_output(L):
                print(f"Found output: \"{self.folder}/{self.name}/out_{L}x{L}.txt")

        for L in self.lengths:
            if self.has_output(L):
                results[L] = RustMonteCarloData(f"{self.folder}/{self.name}/out_{L}x{L}.txt")
                print(f">> Elapsed time for size {L}x{L}: {results[L].elapsed_time/60:.3}min ({results[L].elapsed_time}s)")
                if self.temperatures is None:
                    self.temperatures = results[L].temperatures
                    temp_set = True
        if temp_set:
            print(">> Note: temperatures set from output file.")
        return results