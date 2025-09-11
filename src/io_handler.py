import numpy as np
import subprocess

def perform_rust_computation(parameter_file: str) -> None:
    command = f"cargo run --release -- {parameter_file}"
    
    print(f"Command: \"{command}\"")   
    with subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as stream:
        for line in stream.stdout:
            print(line, end='') 
    
        stdout, stderr = stream.communicate()
        if stream.returncode != 0:
            print(stdout)
            print(stderr)
        else:
            print("=> execution successful.")
    print("")

def array_to_str(array: np.ndarray, rounding: int) -> list[str]:
    temps_str = []
    for val in array:
        temps_str.append(str(round(val,rounding)))
    return ", ".join(temps_str)

def write_parameter_file(param_file: str, out_file: str, rows: np.uint, cols: np.uint, therm_steps: np.uint, measure_steps: np.uint, temperatures: np.uint, rounding = 3)->None:
    # if sys.platform == "win32":
    #     param_file = param_file.replace("/","\\")
    with open(param_file, "w") as f:
        f.write(f"rows: {rows}\n")
        f.write(f"cols: {cols}\n")
        f.write(f"therm_steps: {therm_steps}\n")
        f.write(f"measure_steps: {measure_steps}\n")
        f.write(f"temperatures: {array_to_str(temperatures, rounding)}\n")
        f.write(f"outputfile: {out_file}\n")
        
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
        
        self.all_lists_to_array()

        
if __name__ == "__main__":
    
    data = RustMonteCarloData("no_file_test")