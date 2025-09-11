use std::{usize};
use rand::{rngs, Rng, SeedableRng};
use rand::rngs::SmallRng;
use std::env;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use swendsen_wang::{MonteCarloResults, arange};
use rayon::{self};

use swendsen_wang::{ClusterLabels, IsingArray2D,EquivalenceClass};


fn should_add_to_cluster(rng: &mut SmallRng, J_int: f64, temp: f64) -> bool
{
    let p = 1f64 - (-2_f64*J_int/temp).exp();
    rng.random::<f64>() < p
}

struct SwendsenWangAlgorithm
{
    labels: ClusterLabels,
    eq_classes: EquivalenceClass,
    cluster_flip_probabilities: Vec<f32>
    // cluster_flip_probabilities: Vec<Option<bool>>
}
impl SwendsenWangAlgorithm
{
    fn new(rows: usize, cols: usize) -> Self
    {
        let labels = ClusterLabels::new(rows, cols);
        let eq_classes = EquivalenceClass::new(rows*cols);
        let cluster_flip_probabilities: Vec<f32>  = vec![Default::default(); rows*cols];
        // let cluster_flip_probabilities  = vec![None; rows*cols];
        Self { labels, eq_classes, cluster_flip_probabilities}
    }
    fn reset_cluster_flip_probabilities(&mut self)
    {
        self.cluster_flip_probabilities.fill(Default::default());
    }
    // fn handle_top_left(&mut self, spins: &IsingArray2D, J_int: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let pos = (0, 0);

    //     let s= spins.at_pos(pos);

    //     let classe = self.eq_classes.create_class();    
    //     self.labels.set(pos, classe);

    //     let spin_product = s*(spins.at(pos.0+1, pos.1) + spins.at(pos.0, pos.1+1));
    //     energy_total += -J_int*(spin_product as f64);
    //     spin_sum += s as f64;
        
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_left_edge(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;

    //     let (Ly, _) = spins.shape();
    //     let x = 0;
    //     for y in 1..Ly-1
    //     {
    //         let pos= (y,x);
    //         let above= (y-1,x) ;

    //         let s           = spins.at_pos(pos);
    //         let above_spin  = spins.at_pos(above); // ex: better handle boundary terms, take it out of the loop
            
    //         if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //         {
    //             let below_label = self.labels.at_pos(above);
    //             self.labels.set(pos, below_label);
    //         }
    //         else
    //         {
    //             let new_class = self.eq_classes.create_class();    
    //             self.labels.set(pos, new_class);
    //         }
    //         let spin_product = s*(spins.at(y+1, x) + spins.at(y, x+1));
    //         energy_total += -J_int*(spin_product as f64);
    //         spin_sum += s as f64;
    //     }


    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_bottom_left(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let (Ly, Lx) = spins.shape();

    //     let pos  = (Ly-1,0);
    //     let top  = (0,0);
    //     let above= (Ly-2,0) ;

    //     let s = spins.at_pos(pos);

    //     let above_spin = spins.at_pos(above); // for ex: better handle boundary terms, take it out of the loop
    //     let top_spin   = spins.at_pos(top);
        
        
    //     if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //     {
    //         let below_label = self.labels.at_pos(above);
    //         let cluster_label = self.eq_classes.find(below_label);
    //         self.labels.set(pos, cluster_label);
    //     }
    //     else
    //     {
    //         let new_class = self.eq_classes.create_class();    
    //         self.labels.set(pos,new_class);
    //     }

    //     // PBC:
    //     if s == top_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         let top_label = self.labels.at_pos(top);
    //         let bottom_label = self.labels.at_pos(pos);
    //         let label = self.eq_classes.union_get_label(top_label, bottom_label);
    //         self.labels.set(pos, label);
    //     }

    //     let spin_product = s*(spins.at_periodic(pos.0+1, pos.1) + spins.at(pos.0, pos.1+1));
    //     energy_total += -J*(spin_product as f64);
    //     spin_sum += s as f64;
    
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_bottom(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let Lx = spins.cols;
    //     let Ly = spins.rows;
        
    //     let y = Ly-1;
    //     for x in 1..Lx-1
    //     {
    //         let pos= (y,x);
    //         let s = spins.at_pos(pos);


    //         let left   = (y,x-1);
    //         let above  = (y-1,x) ;
    //         let top    = (0,x);
    //         let bottom = (Ly-1,x);

    //         let left_spin      = spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //         let abvove_spin    = spins.at_pos(above); // ex: better handle boundary terms, take it out of the loop
    //         let top_spin       = spins.at_pos(top);
            
    //         if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             let left_label = self.labels.at_pos(left);
    //             if s == abvove_spin && should_add_to_cluster(rng, J_int, temp)
    //             {
    //                 let below_label = self.labels.at_pos(above);
    //                 let smallest_label = self.eq_classes.union_get_label(left_label, below_label);
    //                 self.labels.set(pos, smallest_label);
    //             }
    //             else
    //             {
    //                 self.labels.set(pos, left_label);
    //             }
    //         }
    //         else if s == abvove_spin && should_add_to_cluster(rng, J_int, temp) 
    //         {
    //             let below_label = self.labels.at_pos(above);
    //             self.labels.set(pos, below_label);
    //         }
    //         else
    //         {
    //             self.eq_classes.create_class();    
    //             self.labels.set(pos, self.eq_classes.data[0]);
    //         }

    //         // connect bottom cluster & top cluster:
    //         if s == top_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             let top_label = self.labels.at_pos(top);
    //             let bottom_label = self.labels.at_pos(bottom);
    //             let label = self.eq_classes.union_get_label(top_label, bottom_label);
    //             self.labels.set(pos, label);
    //         }

    //         let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
    //         energy_total += -J_int*(spin_product as f64);
    //         spin_sum += s as f64;
    //     }
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }

//          TOP
    /////////////////
    //(0,0) (0,1) ....
    //(1,0) (1,1) ....
    //
    /////////////////
//        BOTTOM


    fn perform_swendsen_wang(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    {
        let mut energy_total = 0_f64;
        let mut spin_sum = 0_f64;
        let (Ly, Lx) = spins.shape();
        
        for y in spins.rows()
        {
            for x in spins.columns()                                 
            {
                let pos= (y,x);
                let s = spins.at_pos(pos);


                let left       = (y,x-1);
                let above      = (y-1,x) ;
                let right_edge = (y, Lx-1);
                let left_edge  = (y,0);
                let top    = (0,x);
                let bottom = (Ly-1,x);

                let left_spin      =  if x > 0 {spins.at_pos(left)} else {0}; // NEED TO REDO TO AVOID IF STATEMENT
                let above_spin     =  if y > 0 {spins.at_pos(above)} else {0}; // for ex: better handle boundary terms, take it out of the loop
                let top_spin       = spins.at_pos(top);
                let left_edge_spin = spins.at_pos(left_edge);
                
                if s == left_spin && should_add_to_cluster(rng, J_int, temp)
                {
                    let left_label = self.labels.at_pos(left);
                    if s == above_spin && should_add_to_cluster(rng, J_int, temp)
                    {
                        let below_label = self.labels.at_pos(above);
                        let smallest_label = self.eq_classes.union_get_label(left_label, below_label);
                        self.labels.set(pos, smallest_label);
                    }
                    else
                    {
                        let cluster_label = self.eq_classes.find(left_label);
                        self.labels.set(pos, cluster_label);
                    }
                }
                else if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
                {
                    let below_label = self.labels.at_pos(above);
                    let cluster_label = self.eq_classes.find(below_label);
                    self.labels.set(pos, cluster_label);
                }
                else
                {
                    let new_class = self.eq_classes.create_class();    
                    self.labels.set(pos,new_class);
                }

                // PBC:
                if pos == bottom && s == top_spin && should_add_to_cluster(rng, J_int, temp)
                {
                    let top_label = self.labels.at_pos(top);
                    let bottom_label = self.labels.at_pos(bottom);
                    let label = self.eq_classes.union_get_label(top_label, bottom_label);
                    self.labels.set(pos, label);
                }
                if pos == right_edge && s == left_edge_spin && should_add_to_cluster(rng, J_int, temp)
                {
                    let right_label = self.labels.at_pos(right_edge);
                    let left_label = self.labels.at_pos(left_edge);
                    let label = self.eq_classes.union_get_label(right_label, left_label);
                    self.labels.set(pos, label);
                }

                let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
                energy_total += -J*(spin_product as f64);
                spin_sum += s as f64;
            }
        }
        spin_sum = spin_sum.abs();
        return (energy_total, spin_sum)
    }

    fn flip_cluster(&mut self, spins: &mut IsingArray2D, rng: &mut SmallRng)
    {
        let p_flip = 0.5_f32;
        for i in spins.rows()
        {
            for j in spins.columns()
            {
                let pos = (i,j);
                let label = self.labels.at_pos(pos);
                let cluster_class = self.eq_classes.find(label)-1; // We start at label = 1

                if self.cluster_flip_probabilities[cluster_class] == 0_f32
                {
                    self.cluster_flip_probabilities[cluster_class] = rng.random();
                }
                if self.cluster_flip_probabilities[cluster_class] < p_flip
                {
                    spins.flip_at(i, j);
                }
            }
        }
    }
    fn reset(&mut self)
    {
        self.reset_cluster_flip_probabilities();
        self.eq_classes.reset();
        self.labels.reset();
    }
}

fn perform_swendsen_wang_monte_carlo(J_int: f64, rows: usize, cols: usize, temperatures: Vec<f64>, therm_steps: usize, measure_steps: usize) -> Vec<MonteCarloResults<f64>>
{
    let mut results = vec![MonteCarloResults::<f64>::default(); temperatures.len()];
    // temperatures.iter().zip(&mut results).for_each(|(&temp, result)|
    (&temperatures, &mut results).into_par_iter().for_each(|(&temp, result)| 
    {
        let mut rng   = rngs::SmallRng::from_os_rng();
        let mut spins = IsingArray2D::new_randomized(&mut rng, rows, cols);
        let mut swendsen_wang = SwendsenWangAlgorithm::new(rows, cols);
        for _ in 0..therm_steps
        {
            swendsen_wang.perform_swendsen_wang(&mut spins, &mut rng, J_int, temp);
            swendsen_wang.flip_cluster(&mut spins, &mut rng);
            swendsen_wang.reset();
        }

        let mut energy_acc     = 0_f64;
        let mut energy_sqr_acc = 0_f64;
        let mut spin_sum_acc   = 0_f64;
        let mut spin_sqr_acc   = 0_f64;
        for _ in 0..measure_steps
        {
            let (energy, spin_sum) = swendsen_wang.perform_swendsen_wang(&mut spins, &mut rng, J_int, temp);
            swendsen_wang.flip_cluster(&mut spins, &mut rng);
            swendsen_wang.reset();

            energy_acc     += energy;
            energy_sqr_acc += energy*energy;
            spin_sum_acc   += spin_sum;
            spin_sqr_acc   += spin_sum*spin_sum;
        }

        result.spins_sum_avg = spin_sum_acc/(measure_steps as f64);
        result.spins_sqr_avg = spin_sqr_acc/(measure_steps as f64);
        result.energy_avg    = energy_acc/(measure_steps as f64);
        result.energy_sqr_avg = energy_sqr_acc/(measure_steps as f64);

        // dbg!(spins);
    });

    results
}

const J: f64 = 1_f64;
fn main() 
{
    let args = env::args().collect::<Vec<_>>();
    let usage = "cargo run ROWS COLS --release";
    if args.len() < 5
    {
        print!("{usage}");
        std::process::exit(1);
    }
    let rows: usize = args[1].parse().expect(usage);
    let cols: usize = args[2].parse().expect(usage);
    let therm_steps: usize  = args[3].parse().expect(usage);
    let measure_steps: usize = args[4].parse().expect(usage);


    let temp_start = 1_f64;
    let temp_stop = 3_f64;
    let step = 0.05;

    println!("Launching for {rows}x{cols} with therm steps {therm_steps} & measure_steps: {measure_steps}");
    let temperatures = arange(temp_start, temp_stop, step).expect("step has same sign as (stop-start)");

    let time =     std::time::SystemTime::now();
    let results: Vec<MonteCarloResults<f64>> = perform_swendsen_wang_monte_carlo(J, rows, cols, temperatures.clone(), therm_steps, measure_steps);
    let elapsed_time = time.elapsed().unwrap();
    println!("Time taken: {}s", elapsed_time.as_secs());
    
    let &last_time = temperatures.last().unwrap();
    let file_name = format!("results/{rows}x{cols}_from_{temp_start:.3}_to_{last_time:.3}_with_step_{step:.3}.txt");
    MonteCarloResults::write_to_file(&file_name, &temperatures, &results, rows*cols, elapsed_time).unwrap_or_else(|err|
    {
        print!("Error: {err}");
        std::process::exit(1);
    });
    println!("File saved as {file_name}");

}   
