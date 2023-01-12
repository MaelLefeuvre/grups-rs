pub mod genotype_reader;

mod snpreader;
pub use snpreader::{SNPReader, SNPReaderError};

mod panel_reader;
pub use panel_reader::PanelReader;

mod sampletag;
pub use sampletag::SampleTag;