mod pedparam;
use pedparam::PedParam;

mod param_rate_generator;
pub use param_rate_generator::ParamRateGenerator;
#[allow(clippy::module_inception)]
mod pedigree_params;
pub use pedigree_params::PedigreeParams;