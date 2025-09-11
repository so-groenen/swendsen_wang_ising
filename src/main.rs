#![allow(non_snake_case)]

use std::{usize};
use rand::{rngs, Rng, SeedableRng};
use rand::rngs::SmallRng;
use std::env;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rayon::{self};

use swendsen_wang::{ClusterLabels, IsingArray2D,EquivalenceClass, MonteCarloResults};
use swendsen_wang::parameter_reader;

fn should_add_to_cluster(rng: &mut SmallRng, J_int: f64, temp: f64) -> bool
{
    let p: f64 = 1f64 - (-2_f64*J_int/temp).exp();
    rng.random::<f64>() < p
}
const J: f64 = 1_f64;
trait Direction 
{
    fn left(self) -> Self;
    fn right(self)-> Self;
    fn above(self) -> Self;
    fn below(self) -> Self;
}

impl Direction for (i32, i32)
{
    #[inline(always)]
    fn left(self) -> Self 
    {
        (self.0, self.1 - 1)
    }    
    #[inline(always)]
    fn right(self) -> Self 
    {
        (self.0, self.1 + 1)
    }    
    #[inline(always)]
    fn above(self) -> Self 
    {
        (self.0-1, self.1)
    }    
    #[inline(always)]
    fn below(self) -> Self
    {
        (self.0+1, self.1)    
    }
}
struct SwendsenWangAlgorithm
{
    labels: ClusterLabels,
    eq_classes: EquivalenceClass,
    cluster_flip_probabilities: Vec<f32>
}
impl SwendsenWangAlgorithm
{
    fn new(rows: usize, cols: usize) -> Self
    {
        let labels = ClusterLabels::new(rows, cols);
        let eq_classes = EquivalenceClass::new(rows*cols);
        let cluster_flip_probabilities: Vec<f32>  = vec![Default::default(); rows*cols];
        Self { labels, eq_classes, cluster_flip_probabilities}
    }
    #[inline(always)]
    fn reset_cluster_flip_probabilities(&mut self)
    {
        self.cluster_flip_probabilities.fill(Default::default());
    }
    // fn merge_clusters_single(&mut self, pos: (i32, i32), other: (i32, i32))
    // {
    //     let other_label = self.labels.at_pos(other);
    //     let cluster_label = self.eq_classes.find(other_label);
    //     self.labels.set(pos, cluster_label);
    // }
    #[inline(always)]
    fn merge_clusters_above(&mut self, pos: (i32, i32), above: (i32, i32))
    {
        let above_label = self.labels.at_pos(above);
        let cluster_label = self.eq_classes.find(above_label);
        self.labels.set(pos, cluster_label);
    }
    #[inline(always)]
    fn merge_clusters_left(&mut self, pos: (i32, i32), left_label: usize)
    {
        let cluster_label = self.eq_classes.find(left_label);
        self.labels.set(pos, cluster_label);
    }
    #[inline(always)]
    fn merge_clusters_above_and_left(&mut self, pos: (i32, i32), above: (i32, i32), left_label: usize)
    {
        let above_label = self.labels.at_pos(above); // for ex: above & to the left
        let label = self.eq_classes.union_get_label(above_label, left_label);
        self.labels.set(pos, label);
    }
    #[inline(always)]
    fn make_new_cluster(&mut self, pos: (i32, i32))
    {
        let new_class = self.eq_classes.create_class();    // Each spin is their own class, we use the first element to store the number of classes
        self.labels.set(pos,new_class);        
    }
    #[inline(always)]
    fn merge_clusters_pbc(&mut self, pos: (i32, i32), other_side: (i32, i32))
    {
        let other_side_label = self.labels.at_pos(other_side);
        let current_label = self.labels.at_pos(pos);
        let label = self.eq_classes.union_get_label(current_label, other_side_label);
        self.labels.set(pos, label);
    }
    fn perform_swendsen_wang(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
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

                let left       = pos.left();//(y,x-1);       // I use notations (like in graphics programming) where (i,j)=(y,x)=(0,0) is on the top left
                let above      = pos.above();//(y-1,x) ;      // and the position on the screen decreases as one increase y
                let right_edge = (y, Lx-1);
                let left_edge  = (y,0);
                let top        = (0,x);
                let bottom     = (Ly-1,x);
                 
                let left_spin      =  if x > 0 {spins.at_pos(left)} else {0}; // NEED TO REDO TO AVOID IF STATEMENT
                let above_spin     =  if y > 0 {spins.at_pos(above)} else {0}; // for ex: better handle boundary terms, take it out of the loop
                let top_spin       = spins.at_pos(top);
                let left_edge_spin = spins.at_pos(left_edge);
                
                if s == left_spin && rng.random_bool(proba_add) //should_add_to_cluster(rng, J_int, temp)
                {
                    let left_label = self.labels.at_pos(left);
                    if s == above_spin && rng.random_bool(proba_add)
                    {
                        self.merge_clusters_above_and_left(pos, above, left_label);
                    }
                    else
                    {
                        self.merge_clusters_left(pos, left_label);
                    }
                }
                else if s == above_spin && rng.random_bool(proba_add) 
                {
                    self.merge_clusters_above(pos, above);
                }
                else
                {
                   self.make_new_cluster(pos);
                }

                // Connect bottom & top as well as right edge to left edge:
                if pos == bottom && s == top_spin && rng.random_bool(proba_add)
                {
                    self.merge_clusters_pbc(pos, top);
                }
                if pos == right_edge && s == left_edge_spin && rng.random_bool(proba_add)
                {
                    self.merge_clusters_pbc(pos, left_edge);
                }

                let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
                energy_total += -J*(spin_product as f64);
                spin_sum += s as f64;
            }
        }
        spin_sum = spin_sum.abs();
        return (energy_total, spin_sum)
    }
    // fn handle_top_left(&mut self, spins: &IsingArray2D, J_int: f64) -> (f64, f64)
    // {
    //     let pos = (0, 0);

    //     let s= spins.at_pos(pos);
    //     self.make_new_cluster(pos);

    //     let spin_product = s*(spins.at(pos.0+1, pos.1) + spins.at(pos.0, pos.1+1));
    //     let energy_total= -J_int*(spin_product as f64);
    //     let spin_sum = s as f64;
        
    //     return (energy_total, spin_sum.abs())
    // }
    // fn handle_left_edge(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;

    //     let (Ly, _) = spins.shape();
    //     let x = 0;
    //     for y in 1..Ly-1
    //     {
    //         let pos  = (y,x);
    //         let above= (y-1,x);

    //         let s           = spins.at_pos(pos);
    //         let above_spin  = spins.at_pos(above); // ex: better handle boundary terms, take it out of the loop
            
    //         if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //         {
    //             self.merge_clusters_above(pos, above);
    //         }
    //         else
    //         {
    //             self.make_new_cluster(pos);
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
    //     let (Ly, _) = spins.shape();

    //     let pos  = (Ly-1,0);
    //     let top  = (0,0);
    //     let above= (Ly-2,0);

    //     let s = spins.at_pos(pos);

    //     let above_spin = spins.at_pos(above);  
    //     let top_spin   = spins.at_pos(top);
        
        
    //     if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //     {
    //         self.merge_clusters_above(pos, above);
    //     }
    //     else
    //     {
    //         self.make_new_cluster(pos);
    //     }

    //     // PBC:
    //     if s == top_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         self.merge_clusters_pbc(pos, top);
    //     }

    //     let spin_product = s*(spins.at_periodic(pos.0+1, pos.1) + spins.at(pos.0, pos.1+1));
    //     let energy_total = -J_int*(spin_product as f64);
    //     let spin_sum  = (s as f64).abs();
    
    //     return (energy_total, spin_sum)
    // }
    // fn handle_bottom(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let (Ly, Lx) = spins.shape();
        
    //     let y = Ly-1;

    //     for x in 1..Lx-1                                 
    //     {
    //         let pos= (y,x);
    //         let top =(0,x);
    //         let s = spins.at_pos(pos);


    //         let left       = (y,x-1);
    //         let above      = (y-1,x) ;

    //         let left_spin   = spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //         let above_spin  = spins.at_pos(above); // for ex: better handle boundary terms, take it out of the loop
    //         let top_spin    = spins.at_pos(top); // for ex: better handle boundary terms, take it out of the loop

    //         if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             let left_label = self.labels.at_pos(left);
    //             if s == above_spin && should_add_to_cluster(rng, J_int, temp)
    //             {
    //                 self.merge_clusters_above_and_left(pos, above, left_label);
    //             }
    //             else
    //             {
    //                 self.merge_clusters_left(pos, left_label);
    //             }
    //         }
    //         else if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //         {
    //             self.merge_clusters_above(pos, above);
    //         }
    //         else
    //         {
    //             self.make_new_cluster(pos);
    //         }

    //         // PBC:
    //         if s == top_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             self.merge_clusters_pbc(pos, top);
    //         }

    //         let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
    //         energy_total += -J_int*(spin_product as f64);
    //         spin_sum += s as f64;
    //     }
        
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_bottom_right(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let (Ly, Lx) = spins.shape();
        
    //     let (y,x) = (Ly-1, Lx-1);
    //     let pos= (y,x);
    //     let s = spins.at_pos(pos);


    //     let left       = (y,x-1);
    //     let above      = (y-1,x) ;
    //     let left_edge  = (y,0);
    //     let top        = (0,x);

    //     let left_spin      =  spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //     let above_spin     =  spins.at_pos(above); // for ex: better handle boundary terms, take it out of the loop
    //     let top_spin       = spins.at_pos(top);
    //     let left_edge_spin = spins.at_pos(left_edge);
        
    //     if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         let left_label = self.labels.at_pos(left);
    //         if s == above_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             self.merge_clusters_above_and_left(pos, above, left_label);
    //         }
    //         else
    //         {
    //             self.merge_clusters_left(pos, left_label);
    //         }
    //     }
    //     else if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //     {
    //         self.merge_clusters_above(pos, above);
    //     }
    //     else
    //     {
    //         self.make_new_cluster(pos);
    //     }

    //     // PBC:
    //     if s == top_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         self.merge_clusters_pbc(pos, top);
    //     }
    //     if s == left_edge_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         self.merge_clusters_pbc(pos, left_edge);
    //     }

    //     let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
    //     energy_total += -J_int*(spin_product as f64);
    //     spin_sum += s as f64;
    
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_right_edge(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let (Ly, Lx) = spins.shape();
    //     let x= Lx-1;
    //     for y in 1..Ly-1
    //     {
    //         let pos= (y,x);
    //         let s = spins.at_pos(pos);

    //         let left       = (y,x-1);
    //         let above      = (y-1,x) ;
    //         let left_edge  = (y,0);

    //         let left_spin      =  spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //         let above_spin     =  spins.at_pos(above); // for ex: better handle boundary terms, take it out of the loop
    //         let left_edge_spin = spins.at_pos(left_edge);
            
    //         if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             let left_label = self.labels.at_pos(left);
    //             if s == above_spin && should_add_to_cluster(rng, J_int, temp)
    //             {
    //                 self.merge_clusters_above_and_left(pos, above, left_label);
    //             }
    //             else
    //             {
    //                 self.merge_clusters_left(pos, left_label);
    //             }
    //         }
    //         else if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //         {
    //             self.merge_clusters_above(pos, above);
    //         }
    //         else
    //         {
    //             self.make_new_cluster(pos);
    //         }

    //         // PBC:
    //         if s == left_edge_spin && should_add_to_cluster(rng, J_int, temp)
    //         {
    //             self.merge_clusters_pbc(pos, left_edge);
    //         }

    //         let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
    //         energy_total += -J_int*(spin_product as f64);
    //         spin_sum += s as f64;
    //     }
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }
    // fn handle_top_right(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {

    //     let (Ly, Lx) = spins.shape();
    //     let x = Lx-1;
    //     let y = Ly-1;

    //     let pos= (y,x);
    //     let s = spins.at_pos(pos);

    //     let left       = (y,x-1);
    //     let left_edge  = (y,0);

    //     let left_spin      =  spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //     let left_edge_spin = spins.at_pos(left_edge);
        
    //     if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         let left_label = self.labels.at_pos(left);
    //         self.merge_clusters_left(pos, left_label);
    //     }
    //     else
    //     {
    //         self.make_new_cluster(pos);
    //     }

    //     // PBC:
    //     if s == left_edge_spin && should_add_to_cluster(rng, J_int, temp)
    //     {
    //         self.merge_clusters_pbc(pos, left_edge);
    //     }

    //     let spin_product = s*(spins.at_periodic(y+1, x) + spins.at_periodic(y, x+1));
    //     let energy_total = -J_int*(spin_product as f64);
    //     let spin_sum = (s as f64).abs();
 
    //      return (energy_total, spin_sum)
    // }
    // fn handle_bulk(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let mut energy_total = 0_f64;
    //     let mut spin_sum = 0_f64;
    //     let (Ly, Lx) = spins.shape();
        
    //     for y in 1..Ly-1
    //     {
    //         for x in 1..Lx-1                             
    //         {
    //             let pos= (y,x);
    //             let s = spins.at_pos(pos);

    //             let left       = (y,x-1);
    //             let above      = (y-1,x) ;

    //             let left_spin      = spins.at_pos(left); // NEED TO REDO TO AVOID IF STATEMENT
    //             let above_spin     =  if y > 0 {spins.at_pos(above)} else {0}; // for ex: better handle boundary terms, take it out of the loop
                
    //             if s == left_spin && should_add_to_cluster(rng, J_int, temp)
    //             {
    //                 let left_label = self.labels.at_pos(left);
    //                 if s == above_spin && should_add_to_cluster(rng, J_int, temp)
    //                 {
    //                     self.merge_clusters_above_and_left(pos, above, left_label);
    //                 }
    //                 else
    //                 {
    //                     self.merge_clusters_left(pos, left_label);
    //                 }
    //             }
    //             else if s == above_spin && should_add_to_cluster(rng, J_int, temp) 
    //             {
    //                 self.merge_clusters_above(pos, above);
    //             }
    //             else
    //             {
    //                self.make_new_cluster(pos);
    //             }

    //             let spin_product = s*(spins.at(y+1, x) + spins.at(y, x+1));
    //             energy_total += -J_int*(spin_product as f64);
    //             spin_sum += s as f64;
    //         }
    //     }
    //     spin_sum = spin_sum.abs();
    //     return (energy_total, spin_sum)
    // }

    // fn hoshen_koppelman(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
    // {
    //     let (dE1, dS1) = self.handle_top_left(spins, J_int);
    //     let (dE2, dS2) =self.handle_left_edge(spins, rng, J_int, temp);
    //     let (dE3, dS3) =self.handle_bottom_left(spins, rng, J_int, temp);
    //     let (dE4, dS4) =self.handle_bottom(spins, rng, J_int, temp);
    //     let (dE5, dS5) =self.handle_bottom_right(spins, rng, J_int, temp);
    //     let (dE6, dS6) =self.handle_right_edge(spins, rng, J_int, temp);
    //     let (dE7, dS7) =self.handle_top_right(spins, rng, J_int, temp);
    //     let (dE8, dS8) =self.handle_bulk(spins, rng, J_int, temp);
    //     return (dE1+dE2+dE3+dE4+dE5+dE6+dE7+dE8, dS1+dS2+dS3+dS4+dS5+dS6+dS7+dS8 );
    // }

    
// //          TOP
    /////////////////
    //(0,0) (0,1) ....
    //(1,0) (1,1) ....
    //
    /////////////////
//        BOTTOM


    fn perform_swendsen_wang_explicit(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, J_int: f64, temp: f64) -> (f64, f64)
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


                // let left       = pos.left();
                // let above      = pos.above();
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
                    let above_label = self.labels.at_pos(above);
                    let cluster_label = self.eq_classes.find(above_label);
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
                energy_total += -J_int*(spin_product as f64);
                spin_sum += s as f64;
            }
        }
        spin_sum = spin_sum.abs();
        return (energy_total, spin_sum)
    }

    fn perform_swendsen_wang_explicit2(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
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


                // let left       = pos.left();
                // let above      = pos.above();
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
                
                if s == left_spin && rng.random_bool(proba_add)
                {
                    let left_label = self.labels.at_pos(left);
                    if s == above_spin && rng.random_bool(proba_add)
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
                else if s == above_spin && rng.random_bool(proba_add)//rng.random::<f32>() < proba_add
                {
                    let above_label = self.labels.at_pos(above);
                    let cluster_label = self.eq_classes.find(above_label);
                    self.labels.set(pos, cluster_label);
                }
                else
                {
                    let new_class = self.eq_classes.create_class();    
                    self.labels.set(pos,new_class);
                }

                // PBC:
                if pos == bottom && s == top_spin && rng.random_bool(proba_add)
                {
                    let top_label = self.labels.at_pos(top);
                    let bottom_label = self.labels.at_pos(bottom);
                    let label = self.eq_classes.union_get_label(top_label, bottom_label);
                    self.labels.set(pos, label);
                }
                if pos == right_edge && s == left_edge_spin && rng.random_bool(proba_add)
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

fn perform_swendsen_wang_monte_carlo(rows: usize, cols: usize, temperatures: Vec<f64>, therm_steps: usize, measure_steps: usize) -> Vec<MonteCarloResults<f64>>
{
    let mut results = vec![MonteCarloResults::<f64>::default(); temperatures.len()];
    // temperatures.iter().zip(&mut results).for_each(|(&temp, result)|
    (&temperatures, &mut results).into_par_iter().for_each(|(&temp, result)| 
    {
        let mut rng   = rngs::SmallRng::from_os_rng();
        let mut spins = IsingArray2D::new_randomized(&mut rng, rows, cols);
        let mut swendsen_wang = SwendsenWangAlgorithm::new(rows, cols);

        let proba_add  = 1f64 - (-2_f64*J/temp).exp();

        for _ in 0..therm_steps
        {
            // swendsen_wang.perform_swendsen_wang_explicit(&mut spins, &mut rng, J_int, temp);
            swendsen_wang.perform_swendsen_wang(&mut spins, &mut rng, proba_add);
            swendsen_wang.flip_cluster(&mut spins, &mut rng);
            swendsen_wang.reset();
        }

        let mut energy_acc     = 0_f64;
        let mut energy_sqr_acc = 0_f64;
        let mut spin_sum_acc   = 0_f64;
        let mut spin_sqr_acc   = 0_f64;
        for _ in 0..measure_steps
        {
            // let (energy, spin_sum) = swendsen_wang.perform_swendsen_wang_explicit(&mut spins, &mut rng, J_int, temp);
            let (energy, spin_sum) = swendsen_wang.perform_swendsen_wang(&mut spins, &mut rng, proba_add);
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


const PARAMETERS: [&'static str; 6] = [
    "rows",
    "cols", 
    "therm_steps", 
    "measure_steps",
    "outputfile",
    "temperatures"
];


// const J: f64 = 1_f64;
fn main() 
{
    let args = env::args().collect::<Vec<_>>();
    let params = parameter_reader::read_parameters(&args, &PARAMETERS).unwrap_or_else(|e|
    {
        println!("Failed to read parameters: {e}");
        std::process::exit(1);
    });
    let cols: usize          = params["cols"].parse().expect("Could not parse row.");
    let rows: usize          = params["rows"].parse().expect("could not parse cols.");
    let therm_steps: usize   = params["therm_steps"].parse().expect("could not parse therm_steps.");
    let outputfile: String   = params["outputfile"].parse().expect("could not parse outputfile.");
    let measure_steps: usize = params["measure_steps"].parse().expect("could not parse temperatures.");
    let temperatures: Vec<f64> = params["temperatures"].split(",").map(|t| t.trim().parse().expect("failed parse f64") ).collect();



    println!("Launching Swendsen-Wang simulation for N:{rows}x{cols} with therm steps {therm_steps} & measure_steps: {measure_steps}");
    let &temp_last  = temperatures.last().unwrap();
    let &temp_first = temperatures.first().unwrap();
    let temp_len  = temperatures.len();
    println!("Using: {temp_len} temperatures values from {temp_first} to {temp_last}");


    let time =     std::time::SystemTime::now();
    let results: Vec<MonteCarloResults<f64>> = perform_swendsen_wang_monte_carlo(rows, cols, temperatures.clone(), therm_steps, measure_steps);
    let elapsed_time = time.elapsed().unwrap();
    
    println!("Time taken: {}s", elapsed_time.as_secs());
    

    MonteCarloResults::write_to_file(&outputfile, &temperatures, &results, rows*cols, elapsed_time).unwrap_or_else(|err|
    {
        print!("Error: {err}");
        std::process::exit(1);
    });
    println!("File saved as {outputfile}");

}   
