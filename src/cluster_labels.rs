pub struct ClusterLabels
{
    data: Vec<usize>,
    rows: usize,
    cols: usize,
}

impl ClusterLabels
{
    pub fn new(rows: usize, cols: usize) -> Self
    {
        let default_value = rows*cols; // we could also use 0, but then we would need starting to count classes at 1
        let data = vec![default_value; rows*cols];
        Self {data, rows, cols}
    }
    #[inline]
    pub fn at_pos(&self, pos: (i32, i32)) -> usize
    {
        let i = pos.0 as usize;
        let j = pos.1 as usize;
        self.data[i*self.cols + j]
    }
    #[inline]
    pub fn set(&mut self, pos: (i32, i32), label: usize)
    {
        let i = pos.0 as usize;
        let j = pos.1 as usize;
        self.data[i*self.cols + j] = label;
    }
    // pub fn at_mut(&mut self, i: usize, j: usize) -> &mut usize
    // {
    //     &mut self.data[i*self.cols + j]
    // }
    #[inline]
    pub fn reset(&mut self)
    {
        self.data.fill(Default::default());
    }
}