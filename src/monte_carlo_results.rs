use num_traits::Float;
use std::f64::consts::PI;
use std::io::Write;
use std::iter::zip;

#[derive(Debug, Default, Clone, Copy)]
pub struct MonteCarloResults<T> where T: Float
{
    pub spins_sum_avg: T,   
    pub spins_sqr_avg: T,  
    pub energy_avg: T,       
    pub energy_sqr_avg: T,  
    pub struct_fact_q0: T,
    pub struct_fact_qx: T,      
}


impl<T> MonteCarloResults<T> where T: Float + std::fmt::Display
{
    pub fn write_to_file(file_name: &String, temperatures: &[T], results: &[MonteCarloResults<T>], rows: usize, cols: usize, elapsed_time: std::time::Duration ) -> std::io::Result<()>
    {
        if temperatures.len() != results.len()
        {
            return Err(std::io::Error::other("Results length should match temperature length"));
        }
        if temperatures.iter().any(|t| t.is_zero() || t.is_nan() || t.is_sign_negative())
        {
            return Err(std::io::Error::other("Temperatures should be positive"));
        }
        

        let mut file= std::fs::File::create(file_name)?;
        writeln!(&mut file, "temp, energy_density, magnetisation, specific_heat, susceptibility, correlation length, elapsed_time: {}", elapsed_time.as_secs())?;

        let qx        = 2_f64 * PI / cols as f64;
        let num_spins = T::from(rows*cols).unwrap();

        for  (&temp, res) in zip(temperatures, results)
        {
            let specific_heat   = (res.energy_sqr_avg - res.energy_avg.powi(2) ) / (temp.powi(2) * num_spins);
            let energy_density  = res.energy_avg / num_spins;
            let magnetisation   = res.spins_sum_avg / num_spins;
            let susceptibility  = (res.spins_sqr_avg - res.spins_sum_avg.powi(2)) / (temp * num_spins);
            let mut corr_length = T::zero();
            if !res.struct_fact_q0.is_zero() && res.struct_fact_qx.is_zero()
            {
                let sf_ratio  = (res.struct_fact_q0/res.struct_fact_qx - T::one()).abs().sqrt();
                corr_length   = sf_ratio / T::from(qx).unwrap();
            }
            writeln!(&mut file, "{temp}, {energy_density}, {magnetisation}, {specific_heat}, {susceptibility}, {corr_length}")?;
        }
    
        Ok(())
    }
}

 

pub fn arange<T>(start: T, stop: T, step: T) -> Result<Vec<T>, &'static str> where T: Float
{
    if step == T::zero()
    {
        return Err("Arange Error: Step must be non zero")
    }

    let direction_sign = (stop - start).signum();
    if step.signum() != direction_sign
    {   
        return Err("Arange Error: if stop > (<) start then step must be positive (negative).")
    }   
    let num_of_values = ((stop - start).abs() / step).round().to_usize().unwrap();
    let my_arange  = (0..num_of_values).map(|val| start + step * T::from(val).unwrap() ).collect::<Vec<_>>();
    Ok(my_arange)
}

