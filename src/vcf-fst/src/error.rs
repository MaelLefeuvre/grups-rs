use thiserror::Error;

#[derive(Error, Debug)]
pub enum GenomeFstError {
    #[error("Failed to open target output file")]
    CreateFile(#[source] std::io::Error),

    #[error("Failed to instantiate FST Set Builder")]
    CreateSetBuilder(#[source] fst::Error),

    #[error("The provided file or directory cannot be parsed into a string")]
    DisplayPath,

    #[error("Failed to instantiate threadpool")]
    BuildThreadPool(#[source] rayon::ThreadPoolBuildError),

    //#[error("Failed to instantiate the VCFIndexer")]
    //BuildVcfIndexer,

    //#[error("Failed to generate FSTSet")]
    //BuildFSTSet(#[source] anyhow::Error)


}