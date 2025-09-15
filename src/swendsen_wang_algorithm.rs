#![allow(non_snake_case)]

mod cluster_labels;
mod ising_array_2d;
mod equivalence_class;

use num::complex::Complex64;
use std::f64::consts::PI;
use rand::rngs::SmallRng;
use rand::Rng;

pub use cluster_labels::ClusterLabels;
pub use equivalence_class::EquivalenceClass;
pub use ising_array_2d::IsingArray2D;


pub const J: f64 = 1.0;


#[allow(dead_code)]
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

 

pub struct SwendsenWangAlgorithm
{
    labels: ClusterLabels,
    eq_classes: EquivalenceClass,
    cluster_flip_probabilities: Vec<f32>,
    pub take_fourier_transform: bool,
    fourier_kernels: Vec<Complex64> // [e^{iqx} for x in [0..Lx]] where q = 2pi/Lx
}
impl SwendsenWangAlgorithm
{
    pub fn new(rows: usize, cols: usize) -> Self
    {
        let labels     = ClusterLabels::new(rows, cols);
        let eq_classes = EquivalenceClass::new(rows*cols);

        let cluster_flip_probabilities: Vec<f32>   = vec![Default::default(); rows*cols];
        let take_fourier_transform                 = false;

        let qx  = 2_f64 *PI / cols as f64;
        let fourier_kernels: Vec<Complex64> = (0..cols).map(|x|  (Complex64::i()* qx * (x as f64)).exp() ).collect();

        Self { labels, eq_classes, cluster_flip_probabilities, take_fourier_transform, fourier_kernels}
    }
    #[inline(always)]
    fn reset_cluster_flip_probabilities(&mut self)
    {
        self.cluster_flip_probabilities.fill(Default::default());
    }
    #[inline(always)]
    fn merge_clusters_above(&mut self, pos: (i32, i32), above: (i32, i32))
    {
        let above_label   = self.labels.at_pos(above);
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
        let label       = self.eq_classes.union_get_label(above_label, left_label);
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
        let current_label    = self.labels.at_pos(pos);
        
        let label = self.eq_classes.union_get_label(current_label, other_side_label);
        self.labels.set(pos, label);
    }
    #[inline(always)]
    fn handle_top(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
    {
        let mut energy_total = 0_f64;
        let mut spin_sum     = 0_f64;
        let (_, Lx)          = spins.shape();
        
        // -------------- Top Left -------------- //
        let top_left      = (0,0);
        let top_left_spin = spins.at_pos(top_left);

        self.make_new_cluster(top_left);

        let spin_product = top_left_spin * spins.get_ngbrs_spin_sum(top_left);
        energy_total    += -J*(spin_product as f64);
        spin_sum        += top_left_spin as f64;

        let y = 0;
        // -------------- Top Left -------------- //
        for x in 1..Lx-1
        {
            let pos  = (y,x);
            let left = pos.left();

            let s         = spins.at_pos(pos);
            let left_spin = spins.at_pos(left); 
            
            if s == left_spin && rng.random_bool(proba_add)  
            {
                let left_label = self.labels.at_pos(left);
                self.merge_clusters_left(pos, left_label);
            }
            else
            {
                self.make_new_cluster(pos);
            }
            let spin_product = s * spins.get_ngbrs_spin_sum(pos);
            energy_total    += -J*(spin_product as f64);
            spin_sum        += s as f64;
        }
        // -----------------------------------------------//
        let top_right   = (0, Lx-1);
        let left        = top_right.left();
        let s           = spins.at_pos(top_right);
        let left_spin   = spins.at_pos(left); 

        if s == left_spin && rng.random_bool(proba_add)  
        {
            let left_label = self.labels.at_pos(left);
            self.merge_clusters_left(top_right, left_label);
        }
        else
        {
            self.make_new_cluster(top_right);
        }

        // connect with top left corner:
        if s == top_left_spin && rng.random_bool(proba_add)
        {
            self.merge_clusters_pbc(top_right, top_left);
        }
        let spin_product = s * spins.get_ngbrs_spin_sum_pbc(top_right);
        energy_total += -J*(spin_product as f64);
        spin_sum     += s as f64;
        
        (energy_total, spin_sum)
    }
    #[inline(always)]
    fn handle_rows(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
    {
        let mut energy_total = 0_f64;
        let mut spin_sum     = 0_f64;
        let (Ly, Lx) = spins.shape();

        for y in 1..Ly-1
        {
            // ---------------------- LEFT EDGE ---------------------------//
            let left_edge      = (y,0);
            let above          = left_edge.above(); 
            let left_edge_spin = spins.at_pos(left_edge);
            let above_spin     = spins.at_pos(above);
            if left_edge_spin == above_spin && rng.random_bool(proba_add) 
            {
                self.merge_clusters_above(left_edge, above);
            }
            else
            {
                self.make_new_cluster(left_edge);
            }
            let spin_product = left_edge_spin * spins.get_ngbrs_spin_sum(left_edge);
            energy_total    += -J*(spin_product as f64);
            spin_sum        += left_edge_spin as f64;
            // ----------------------- BULK -----------------------//
            for x in 1..Lx-1
            {
                let pos   = (y,x);
                let left  = pos.left();
                let above = pos.above(); 
                let s          = spins.at_pos(pos);
                let left_spin  = spins.at_pos(left);  
                let above_spin = spins.at_pos(above);  
                if s == left_spin && rng.random_bool(proba_add)  
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
                let spin_product = s * spins.get_ngbrs_spin_sum(pos);
                energy_total += -J*(spin_product as f64);

                // energy_total += spins.get_energy(J, pos);
                spin_sum     += s as f64;
            }
            // ----------------------- RIGHT EDGE -----------------------//
            let right_edge  = (y, Lx-1);
            let left        = right_edge.left();
            let above       = right_edge.above();
            let s           = spins.at_pos(right_edge);
            let left_spin   = spins.at_pos(left);  
            let above_spin  = spins.at_pos(above);
            if s == left_spin && rng.random_bool(proba_add)  
            {
                let left_label = self.labels.at_pos(left);
                if s == above_spin && rng.random_bool(proba_add)
                {
                    self.merge_clusters_above_and_left(right_edge, above, left_label);
                }
                else
                {
                    self.merge_clusters_left(right_edge, left_label);
                }
            }
            else if s == above_spin && rng.random_bool(proba_add) 
            {
                self.merge_clusters_above(right_edge, above);
            }
            else
            {
                self.make_new_cluster(right_edge);
            }
            // connect with left edge:
            if s == left_edge_spin && rng.random_bool(proba_add)
            {
                self.merge_clusters_pbc(right_edge, left_edge);
            }
            let spin_product = s * spins.get_ngbrs_spin_sum_pbc(right_edge);
            energy_total    += -J*(spin_product as f64);
            spin_sum        += s as f64;
        }
        return (energy_total, spin_sum);
    }
    #[inline(always)]
    fn handle_bottom(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
    {
        let mut energy_total = 0_f64;
        let mut spin_sum     = 0_f64;
        let (Ly, Lx)         = spins.shape();
        
        // ---------------------- Bottom Left ---------------------------//
        let bottom_left      = (Ly-1,0);
        let above            = bottom_left.above(); 
        let top_left         = (0, 0);
        let bottom_left_spin = spins.at_pos(bottom_left);
        let above_spin       = spins.at_pos(above);
        let top_left_spin    = spins.at_pos(top_left);

        if bottom_left_spin == above_spin && rng.random_bool(proba_add) 
        {
            self.merge_clusters_above(bottom_left, above);
        }
        else
        {
            self.make_new_cluster(bottom_left);
        }
        // connect to top left
        if bottom_left_spin == top_left_spin && rng.random_bool(proba_add)
        {
            self.merge_clusters_pbc(bottom_left, top_left);
        }
        let spin_product = bottom_left_spin * spins.get_ngbrs_spin_sum_pbc(bottom_left);
        energy_total    += -J*(spin_product as f64);
        spin_sum        += bottom_left_spin as f64;

        // ---------------------- Bottom ---------------------------//
        for x in 1..Lx-1
        {
            let pos        = (Ly-1, x);
            let top        = (0,    x);
            let left       = pos.left();
            let above      = pos.above();
            let s          = spins.at_pos(pos);
            let left_spin  = spins.at_pos(left); 
            let top_spin   = spins.at_pos(top);
            let above_spin = spins.at_pos(above);

            if s == left_spin && rng.random_bool(proba_add)  
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
            // PBC:
            if s == top_spin && rng.random_bool(proba_add)
            {
                self.merge_clusters_pbc(pos, top);
            }
            let spin_product = s * spins.get_ngbrs_spin_sum_pbc(pos);
            energy_total += -J*(spin_product as f64);
            // energy_total += spins.get_energy_periodic(J, pos);
            spin_sum     += s as f64;
        }

        // ---------------------- Bottom Right ---------------------------//
        let bottom_right   = (Ly-1, Lx-1);
        let top_right      = (0,    Lx-1);
        let left           = bottom_right.left();
        let above          = bottom_right.above();
        let s              = spins.at_pos(bottom_right);
        let left_spin      = spins.at_pos(left); 
        let top_right_spin = spins.at_pos(top_right);
        let above_spin     = spins.at_pos(above);

        if s == left_spin && rng.random_bool(proba_add)  
        {
            let left_label = self.labels.at_pos(left);
            if s == above_spin && rng.random_bool(proba_add)
            {
                self.merge_clusters_above_and_left(bottom_right, above, left_label);
            }
            else
            {
                self.merge_clusters_left(bottom_right, left_label);
            }
        }
        else if s == above_spin && rng.random_bool(proba_add)
        {
            self.merge_clusters_above(bottom_right, above);
        }
        else
        {
            self.make_new_cluster(bottom_right);
        }
        // connect with bottom left corner && top right corner:
        if s == bottom_left_spin && rng.random_bool(proba_add)
        {
            self.merge_clusters_pbc(bottom_right, bottom_left);
        }
        if s == top_right_spin && rng.random_bool(proba_add)
        {
            self.merge_clusters_pbc(bottom_right, top_right);
        }
        let spin_product = s * spins.get_ngbrs_spin_sum_pbc(bottom_right);
        energy_total += -J*(spin_product as f64);
        // energy_total += spins.get_energy_periodic(J, top_right);
        spin_sum     += s as f64;
        
        (energy_total, spin_sum)
    }
    #[inline(always)]
    pub fn perform_swendsen_wang_all(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
    {
        let (dE1, dS1) = self.handle_top(spins, rng, proba_add);
        let (dE2, dS2) = self.handle_rows(spins, rng, proba_add);
        let (dE3, dS3) = self.handle_bottom(spins, rng, proba_add);

        let spin_sum_abs = (dS1+dS2+dS3).abs();
        let total_energy = dE1+dE2+dE3;

        (total_energy, spin_sum_abs)
    }    
    pub fn flip_cluster_and_take_fourier(&mut self, spins: &mut IsingArray2D, rng: &mut SmallRng) -> (f64, Complex64)
    {
        let p_flip      = 0.5_f32;
        let (Ly, Lx)    = spins.shape();
        let factor      = 1_f64 / ((Lx*Ly) as f64).sqrt();
        let mut spin_q0 = f64::default();
        let mut spin_qx = Complex64::default();

        for y in spins.rows()
        {
            for x in spins.columns()
            {
                let pos           = (y,x);
                let label         = self.labels.at_pos(pos);
                let cluster_class = self.eq_classes.find(label)-1; 

                if self.cluster_flip_probabilities[cluster_class] == 0_f32
                {
                    self.cluster_flip_probabilities[cluster_class] = rng.random();
                }
                if self.cluster_flip_probabilities[cluster_class] < p_flip
                {
                    spins.flip_at(y, x);
                }

                if self.take_fourier_transform
                {
                    let s       = spins.at_pos(pos) as f64;
                    let exp_iqx = self.fourier_kernels[x as usize]; //kernels = [e^{i*qx*x} for x in [0..Lx]] where qx = 2pi/Lx
                    spin_q0 += factor * s;
                    spin_qx += factor * s * exp_iqx;
                }
            }
        }
        (spin_q0, spin_qx)
    }
    pub fn reset(&mut self)
    {
        self.reset_cluster_flip_probabilities();
        self.eq_classes.reset();
        self.labels.reset();
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    #[allow(dead_code)] // This version was first adapted from my Python-Numba version. a bit less efficient but more readable!
    fn perform_swendsen_wang_complact_version(&mut self, spins: &IsingArray2D, rng: &mut SmallRng, proba_add: f64) -> (f64, f64)
    {
        let mut energy_total = 0_f64;
        let mut spin_sum     = 0_f64;
        let (Ly, Lx) = spins.shape();
        
        for y in spins.rows()
        {
            for x in spins.columns()                                 
            {
                let pos = (y,x);
                let s   = spins.at_pos(pos);

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
                spin_sum     += s as f64;
            }
        }
        spin_sum = spin_sum.abs();
        return (energy_total, spin_sum)
    }
}
