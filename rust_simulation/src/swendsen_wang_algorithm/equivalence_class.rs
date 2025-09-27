
pub struct EquivalenceClass
{
    data: Vec<usize>,
}
impl EquivalenceClass
{
    pub fn new(num_spins: usize) -> Self
    {
        let data: Vec<usize> = (0..(num_spins+1)).collect();
        Self { data}//, n_classes, }
    }   
    #[inline]
    pub fn reset(&mut self)
    {
        self.data.iter_mut().enumerate().for_each(|(count, value)| *value = count);
    }
    pub fn create_class(&mut self) -> usize
    {       
        self.data[0]         += 1;
        let new_cluster_index = self.data[0];
        
        self.data[new_cluster_index] = new_cluster_index;
        self.data[0]

    }
    pub fn find(&self, mut cluster_result: usize) -> usize 
    {
        while self.data[cluster_result] != cluster_result
        {
            cluster_result = self.data[cluster_result]
        };
        return cluster_result;
    }
    pub fn union_get_label(&mut self, cluster_label1: usize, cluster_label2: usize) -> usize
    {
        let l1    = self.find(cluster_label1);
        let l2    = self.find(cluster_label2);
        let l_max = l1.max(l2);
        let l_min = l1.min(l2);
        if l_min != l_max
        {
            self.data[l_max] = l_min;
        }
        return l_min;
    }
}
