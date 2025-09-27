#![allow(non_snake_case)]
use rand::rngs::SmallRng;
use rand::Rng;

pub struct IsingArray2D
{
    data: Vec<i8>,
    rows: i32,
    cols: i32,
}


trait MonteCarloModulo 
{
    fn modulo(self, other: Self) -> Self;   
}

// Se:: https://rust.godbolt.org/z/x4KMMhn69
impl MonteCarloModulo for i32
{
    #[inline(always)]
    fn modulo(self, other: Self) -> Self
    {
        // (self % other + other) % other
        if self >= other
        {
            return self - other;
        }
        if self < 0
        {
            return self + other;
        }
        self
    }    
}

impl IsingArray2D
{
    pub fn new_polarized(rows: usize, cols: usize) -> Self
    {
        let data = vec![1_i8; rows*cols];
        let rows = rows as i32;
        let cols = cols as i32;

        Self {data, rows, cols}
    }
    pub fn new_randomized(rng: &mut SmallRng, rows: usize, cols: usize) -> Self
    {
        let mut spins = IsingArray2D::new_polarized(rows, cols);
        spins.randomize_spins(rng);
        spins
    }
    #[inline(always)]
    pub fn at_pos(&self, pos: (i32, i32)) -> i8 
    {
        self.at(pos.0, pos.1)
    }
    #[inline(always)]
    pub fn at_pos_periodic(&self, pos: (i32, i32)) -> i8 
    {
        self.at_periodic(pos.0, pos.1)
    }
    #[inline(always)]
    pub fn at(&self, i: i32, j: i32) -> i8
    {
        self.data[(i*self.cols) as usize + j as usize]
    }        
    // #[inline(always)]
    // pub fn get_energy(&self, J_int: f64, pos: (i32, i32)) -> f64
    // {
    //     let s = self.at_pos(pos);
    //     let spin_product = s*(self.at(pos.0 + 1, pos.1) + self.at(pos.0, pos.1+1));
    //     -J_int * (spin_product as f64)
    // }
    #[inline(always)]
    pub fn get_ngbrs_spin_sum(&self, pos: (i32, i32)) -> i8
    {
        self.at(pos.0 + 1, pos.1) + self.at(pos.0, pos.1+1)
    }
    #[inline(always)]
    pub fn get_ngbrs_spin_sum_pbc(&self, pos: (i32, i32)) -> i8
    {
        self.at_periodic(pos.0 + 1, pos.1) + self.at_periodic(pos.0, pos.1+1)
    }
    #[inline(always)]
    pub fn flip_at(&mut self, i: i32, j: i32)
    {
        self.data[(i*self.cols) as usize + j as usize] *= -1;
    }
    #[inline(always)]
    pub fn at_periodic(&self, i: i32, j: i32) -> i8
    {
        // let i = i.rem_euclid(self.rows);
        // let j = j.rem_euclid(self.cols);
        let i = i.modulo(self.rows);
        let j = j.modulo(self.cols);
        debug_assert!(i >= 0, "at periodic: i in [0,rows[");
        debug_assert!(j >= 0, "at periodic: j in [0,cols[");
        self.at(i,j)
    }
    #[inline(always)]
    pub fn shape(&self) -> (i32, i32)
    {
        (self.rows, self.cols)
    }
    #[inline(always)]
    pub fn rows(&self) -> std::ops::Range<i32>
    {
        0..self.rows
    }
    #[inline(always)]
    pub fn columns(&self) -> std::ops::Range<i32>
    {
        0..self.cols
    }
    pub fn randomize_spins(&mut self, rng: &mut SmallRng)
    {
        for i in 0..self.cols
        {
            for j in 0..self.rows
            {
                if rng.random_bool(0.5)
                {
                    self.flip_at(i, j);
                }
            }
        }
    }
}

use std::fmt;
impl fmt::Debug for IsingArray2D
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result 
    {
        let mut result = vec![String::from("\n[")];
        for i in 0..self.rows-1
        {
            for j in 0..self.cols
            {
                result.push(format!(" {} ", self.at(i, j)));
            }
            result.push(String::from("\n"));
        };
        for j in 0..self.cols
        {
            result.push(format!(" {} ", self.at(self.rows-1, j)));
        }
        result.push(String::from("]"));
        let result = result.concat();
        write!(f, "{}", result)
    }
}