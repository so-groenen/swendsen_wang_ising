from __future__ import annotations
import numpy as np

class RustMonteCarloData:
    def __init__(self, file_name):
        self.file_name = file_name
        self.temperatures       = []
        self.energy_density     = []
        self.magnetisation      = []
        self.specific_heat      = []
        self.mag_susceptibility = []
        self.elapsed_time       = -1
        self.observables        = []
        self.correlation_length = []
        self.compute_observables()
        
    def to_nd_array(self, name):
        my_lists = getattr(self, name)
        
        if not isinstance(my_lists, np.ndarray):            
            setattr(self, name, np.asarray(my_lists))

    def all_lists_to_array(self):
        for name in vars(self).keys():
            if name != "elapsed_time":
                self.to_nd_array(name)

    def has_file(self)->bool:
        try:
            with open(self.file_name, "r") as _:
                pass
            return True
        except Exception as e:
            print(f"Could not open file: {e}")
            return False
        
    def compute_observables(self) -> None:
        if not self.has_file():
            return
        with open(self.file_name, "r") as file:
            for (n, lines) in enumerate(file):
                if n == 0:
                    slines = lines.split(':')
                    try:
                        self.observables  = slines[0].split(', ')
                        self.elapsed_time = float(slines[1])
                    except Exception as _:
                        self.observables = lines.split(',')
                        print("No elasped time found.")
                else:
                    slines = lines.split(", ")
                    self.temperatures.append(float(slines[0]))
                    self.energy_density.append(float(slines[1]))
                    self.magnetisation.append(float(slines[2]))
                    self.specific_heat.append(float(slines[3]))
                    self.mag_susceptibility.append(float(slines[4]))
                    self.correlation_length.append(float(slines[5]))

        self.all_lists_to_array()

    @staticmethod
    def have_same_temps(data: dict[int, RustMonteCarloData]):
        for res1 in data.values():
            for res2 in data.values():
                if res1 is not res2 and np.all(res1.temperatures == res2.temperatures):
                    return False
        return False
        
if __name__ == "__main__":
    
    data = RustMonteCarloData("no_file_test")