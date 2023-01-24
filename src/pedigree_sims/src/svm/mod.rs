

use located_error::LocatedError;
use anyhow::Result;
use crate::pedigrees::PedigreeReps;

use libsvm::{SvmInit, SvmTrainer, KernelInit, ModelInit};

mod error;
use error::{SVMError, SvmBuilderError};

mod linfa_svm;
mod smartcore_svm;


pub struct LibSvmBuilder<'a> {
    inner   : SvmTrainer,
    features: Option<Vec<Vec<f64>>>,
    labels  : Option<Vec<f64>>,
    label_order: Option<&'a Vec<&'a String>>
}

impl<'a> Default for LibSvmBuilder<'a> {
    fn default() -> Self {
        // ---- Initialize binary SVM with linear kernel.
        let svm_initializer  = SvmInit {
            kernel                : Some(KernelInit::Linear),
            model                 : Some(ModelInit::CSvc { cost: Some(10.) }),
            label_weights         : None,
            termination_eps       : Some(0.1),
            probability_estimates : Some(false),
            shrinking             : Some(true),
            cache_size            : Some(40)
        };
        let svm_trainer = svm_initializer.build()
            .expect("Failed to initialize libsvm SVMTrainer");
        Self { inner: svm_trainer, features: None, labels: None , label_order: None}
    }
}

impl<'a> LibSvmBuilder<'a> {
    pub fn sims(&mut self, pedigrees: &PedigreeReps) -> &mut Self {
        let mut data: Vec<Vec<f64>> = Vec::with_capacity(pedigrees.len());
        data.extend(
            pedigrees.iter()
            .flat_map(|ped| ped.comparisons.iter()
                .map(|cmp| vec![cmp.get_avg_pwd()])
            )
        );
        self.features = Some(data);
        self
    }

    pub fn label_order(&mut self, label_order: &'a Vec<&String>) -> &mut Self {
        self.label_order = Some(label_order);
        self
    }

    pub fn labels(&mut self, pedigrees: &PedigreeReps, label_treshold: &usize) -> &mut Self {
        let Some(label_order) = self.label_order else {
            panic!()
        };
        let mut targets: Vec<f64> = Vec::with_capacity(pedigrees.len());
        targets.extend(
            pedigrees.iter().flat_map(|ped| 
                ped.comparisons.iter().map(|cmp| {
                    let label_idx = label_order.iter().position(|lab| *lab == &cmp.label )
                        .expect("Invalid Label Found while building SVM labels.");
                    let is_greater_than = label_idx > *label_treshold;
                    (is_greater_than as u32) as f64
                })
            )
        );
        self.labels      = Some(targets);
        self.label_order = Some(label_order);
        self
    }

    pub fn build(&self) -> Result<LibSvm> {
        let loc_msg = "While attempting to build SVM from builder.";
        let Some(ref labels) = self.labels else {
            return Err(SvmBuilderError::MissingLabels).loc(loc_msg)
        };

        let Some(ref features) = self.features else {
            return Err(SvmBuilderError::MissingFeatures).loc(loc_msg)
        };
        
        let svm = self.inner.fit(&features[..], &labels[..])
            .map_err(SvmBuilderError::FitSvm)
            .loc(loc_msg)?;
        
        Ok(LibSvm { inner: svm })
    }
}
pub struct LibSvm {
    inner : libsvm::SvmPredictor,
}

impl LibSvm {
    /// Directly create an SVM from a pedigree vector.
    #[allow(unused)]
    pub fn from_comparisons(pedigrees: &PedigreeReps, label_treshold: &usize, label_order: &Vec<&String>) -> Result<Self> {
    
        let svm = LibSvmBuilder::default()
            .sims(pedigrees)
            .label_order(label_order)
            .labels(pedigrees, label_treshold)
            .build()
            .loc("While attempting to build SVM.")?;

        Ok( Self{inner: svm.inner} )
    }

    #[allow(unused)]
    pub fn predict_probs(&self, avg: f64) -> Result<Vec<(f64, Vec<f64>)>> {
        let x: [f64; 1 ] = [avg];
        let x: &[f64] = x.as_ref();
        self.inner.predict_with_probability(x)
            .map_err(SVMError::PredictProbs)
            .with_loc(|| format!("While attempting to predict value {avg}"))
    }

    pub fn predict(&self, avg: f64) -> Result<f64>{
        let x: [f64; 1 ] = [avg];
        let x: &[f64] = x.as_ref();
        Ok(self.inner.predict(x)
            .map_err(SVMError::PredictBool)
            .with_loc(|| format!("While attempting to predict value {avg}"))?[0]
        )
        
    }
}

