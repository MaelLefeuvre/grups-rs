#[derive(Debug, Clone)]
pub struct SampleTag {
    id : String,
    idx: usize,
}

impl SampleTag {
    pub fn new(id: &str, idx: usize) -> SampleTag {
        SampleTag{id: id.to_string(), idx}
    }

    pub fn id(&self) -> &String {
        &self.id
    }
    pub fn idx(&self) -> &usize {
        &self.idx
    }
}

impl PartialEq for SampleTag {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
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
