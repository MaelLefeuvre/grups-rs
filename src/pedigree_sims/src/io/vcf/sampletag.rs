/// Simple Struct representing an Input Sample id information
/// # Fields
/// - `id` : (String)        - Name of the sample
/// - `idx`: (Option<usize>) - When using a VCFReader, `idx` corresponds to the column field at which the genotype
///                            information  of this sample is located within the vcf file
#[derive(Debug, Clone)]
pub struct SampleTag {
    id : String,
    idx: Option<usize>,
}

impl SampleTag {
    /// Instantiate a new SampleTag struct
    /// # Arguments:
    /// - `id` : raw string slice corresponding to the name of the sample
    /// - `idx`: optional vcf column field index of the sample.
    pub fn new(id: &str, idx: Option<usize>) -> SampleTag {
        SampleTag{id: id.to_string(), idx}
    }

    /// Return the name of the sample.
    pub fn id(&self) -> &String {
        &self.id
    }

    /// Return the vcf column field index of the sample.
    pub fn idx(&self) -> &Option<usize> {
        &self.idx
    }

    /// Set the vcf column field index of the sample.
    pub fn set_idx(&mut self, idx: usize) {
        self.idx = Some(idx);
    }
}

impl PartialEq for SampleTag {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl PartialEq<&str> for SampleTag {
    fn eq(&self, string: &&str) -> bool {
        self.id.as_str() == *string
    }
}

impl Eq for SampleTag {}

impl std::cmp::Ord for SampleTag {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.cmp(&other.id)
    }
}

impl std::cmp::PartialOrd for SampleTag {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::fmt::Display for SampleTag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let idx_str = match self.idx { Some(idx) => format!(" ({})", idx), None => String::from("")};
        write!(f, "{}{}", self.id, idx_str)
    }
}
