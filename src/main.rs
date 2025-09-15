#![allow(non_snake_case)]

use std::{usize};
use num::{Zero};
use rand::{rngs, SeedableRng};
use std::env;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rayon::{self};
use num::complex::ComplexFloat;

use swendsen_wang::monte_carlo_results::MonteCarloResults;
use swendsen_wang::swendsen_wang_algorithm::{SwendsenWangAlgorithm, IsingArray2D};
use swendsen_wang::parameter_reader::ParameterReader;


const MINIMUM_TEMP: f64 = 1E-6;
use swendsen_wang::swendsen_wang_algorithm::J;

 
fn perform_swendsen_wang_monte_carlo(rows: usize, cols: usize, temperatures: Vec<f64>, therm_steps: usize, measure_steps: usize, measure_corr_length: bool) -> Vec<MonteCarloResults<f64>>
{
    let mut results = vec![MonteCarloResults::<f64>::default(); temperatures.len()];

    (&temperatures, &mut results).into_par_iter().for_each(|(&temp, result)| 
    {
        let mut rng           = rngs::SmallRng::from_os_rng();
        let mut spins         = IsingArray2D::new_randomized(&mut rng, rows, cols);
        let mut swendsen_wang = SwendsenWangAlgorithm::new(rows, cols);

        let proba_add = 1f64 - (-2_f64*J/temp).exp();
        
        swendsen_wang.take_fourier_transform = false;
        for _ in 0..therm_steps
        {
            swendsen_wang.perform_swendsen_wang_all(&mut spins, &mut rng, proba_add);
            swendsen_wang.flip_cluster_and_take_fourier(&mut spins, &mut rng);
            swendsen_wang.reset();
        }

        let mut energy_acc     = 0_f64;
        let mut energy_sqr_acc = 0_f64;
        let mut spin_sum_acc   = 0_f64;
        let mut spin_sqr_acc   = 0_f64;
        let mut re_spin_q0_sqr_acc = 0_f64; //  <Re[sigma_q0]²>  
        let mut re_spin_qx_sqr_acc = 0_f64; //  <Re[sigma_qx]²>  
        let mut im_spin_qx_sqr_acc = 0_f64; //  <Im[sigma_qx]²> 

        swendsen_wang.take_fourier_transform = measure_corr_length;
        for _ in 0..measure_steps
        {
            let (energy, spin_sum) = swendsen_wang.perform_swendsen_wang_all(&mut spins, &mut rng, proba_add);
            let (spin_q0, spin_qx) = swendsen_wang.flip_cluster_and_take_fourier(&mut spins, &mut rng);
            swendsen_wang.reset();

            energy_acc     += energy;
            energy_sqr_acc += energy*energy;
            spin_sum_acc   += spin_sum;
            spin_sqr_acc   += spin_sum*spin_sum;

            // Structure factor calculation
            if swendsen_wang.take_fourier_transform
            {
                re_spin_q0_sqr_acc += spin_q0*spin_q0;
                re_spin_qx_sqr_acc += spin_qx.re()*spin_qx.re();
                im_spin_qx_sqr_acc += spin_qx.im()*spin_qx.im();
            }            
        }
        result.struct_fact_q0 = re_spin_q0_sqr_acc/(measure_steps as f64);                          //S(q0) =  <Re[sigma_q0]²>
        result.struct_fact_qx = (re_spin_qx_sqr_acc + im_spin_qx_sqr_acc)/(measure_steps as f64);   //S(qx) =  <Re[sigma_qx]²> + <Im[sigma_qx]²>
        result.spins_sum_avg  = spin_sum_acc/(measure_steps as f64);
        result.spins_sqr_avg  = spin_sqr_acc/(measure_steps as f64);
        result.energy_avg     = energy_acc/(measure_steps as f64);
        result.energy_sqr_avg = energy_sqr_acc/(measure_steps as f64);
    });

    results
}


const PARAMETERS: [&'static str; 7] = [
    "rows",
    "cols", 
    "therm_steps", 
    "measure_steps",
    "outputfile",
    "temperatures",
    "measure_struct_fact"
];

fn main() 
{
    let args   = env::args().collect::<Vec<_>>();
    let reader = ParameterReader::build(&args, &PARAMETERS).unwrap_or_else(|e|
    {
        println!("Failed to create parameter reader: {e}");
        std::process::exit(1);
    });

    let params = reader.parse_parameters(":").unwrap_or_else(|e|
    {
        println!("Failed to read parameters reader: {e}");
        std::process::exit(1);
    });
    
    let cols: usize                = params["cols"].parse().expect("!! Could not parse \"cols\"");
    let rows: usize                = params["rows"].parse().expect("!! Could not parse \"rows\"");
    let therm_steps: usize         = params["therm_steps"].parse().expect("!! Could not parse \"therm_steps\"");
    let outputfile: String         = params["outputfile"].parse().expect("!! Could not parse \"outputfile\"");
    let measure_steps: usize       = params["measure_steps"].parse().expect("!! Could not parse \"measure_steps\"");
    let mut temperatures: Vec<f64> = params["temperatures"].split(",").map(|t| t.trim().parse().expect("!! failed parse temperatures") ).collect();
    let measure_struct_fact: bool  = params["measure_struct_fact"].to_lowercase().parse().expect("!! Could not parse structur factor");

    temperatures
        .iter_mut()
        .filter(|t| t.is_zero() || t.is_nan() || t.is_sign_negative())
        .for_each(|t|
        {
            println!("Setting temperature t={t} => {MINIMUM_TEMP}");
            *t = MINIMUM_TEMP;
        });
    
    let temperatures = temperatures; // -> remove mutability


    println!("Launching Swendsen-Wang simulation for N:{rows}x{cols} with therm steps {therm_steps} & measure_steps: {measure_steps}");
    let &temp_last  = temperatures.last().unwrap();
    let &temp_first = temperatures.first().unwrap();
    let temp_len  = temperatures.len();
    println!("Using: {temp_len} temperatures values from {temp_first} to {temp_last}");
    println!("Measuring correlation length: {measure_struct_fact}");

    let time = std::time::SystemTime::now();
    let results: Vec<MonteCarloResults<f64>> = perform_swendsen_wang_monte_carlo(rows, cols, temperatures.clone(), therm_steps, measure_steps, measure_struct_fact);
    let elapsed_time = time.elapsed().unwrap();
    
    println!("Time taken: {}s", elapsed_time.as_secs());
    

    MonteCarloResults::write_to_file(&outputfile, &temperatures, &results, rows, cols, elapsed_time).unwrap_or_else(|err|
    {
        print!("Could not write to file: {err}");
        std::process::exit(1);
    });
    println!("File saved as {outputfile}");

}   
