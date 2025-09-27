
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


def DEBUG_get_default_monte_carlo_parameters(lengths: list[int]) -> tuple[dict, dict]:
    """Debug params"""   
    thermalisation_steps = dict()
    measurement_steps    = dict()
    for L in lengths:
        thermalisation_steps[L] = int(1e3)
        measurement_steps[L]    = int(1e3)
    return (thermalisation_steps, measurement_steps)


def main():
    print("Hello from simulation-manager!")


if __name__ == "__main__":
    main()
