import numpy as np
import subprocess

def perform_rust_computation(parameter_file: str) -> None:
    command  = f"cargo run --release -- {parameter_file}"
    time     = -1
    
    print(f"Command: \"{command}\"")   
    with subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1, text=True, stderr=subprocess.STDOUT) as stream:
        for line in stream.stdout:
            if line.__contains__("Time taken: "):
                time_str = line.split("Time taken: ")[1]
                time     = int(time_str.removesuffix("s\n"))
                
            print(f" * From Rust: {line}", end='') 

        stdout, _ = stream.communicate()
        if stream.returncode != 0:
            print(stdout)
            # print(stderr)
        else:
            print("=> execution successful.")
    print("")
    return time

def array_to_str(array: np.ndarray, rounding: int) -> list[str]:
    temps_str = []
    for val in array:
        temps_str.append(str(round(val,rounding)))
    return ", ".join(temps_str)

def write_parameter_file(param_file: str, out_file: str, rows: np.uint, cols: np.uint, therm_steps: np.uint, measure_steps: np.uint, temperatures: np.uint, rounding = 3)->None:
    with open(param_file, "w") as f:
        f.write(f"rows: {rows}\n")
        f.write(f"cols: {cols}\n")
        f.write(f"therm_steps: {therm_steps}\n")
        f.write(f"measure_steps: {measure_steps}\n")
        f.write(f"temperatures: {array_to_str(temperatures, rounding)}\n")
        f.write(f"outputfile: {out_file}\n")

if __name__ == "__main__":
    None