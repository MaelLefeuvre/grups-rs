

use located_error::LocatedError;
use anyhow::Result;
use log::{warn, trace};
use crate::pedigrees::PedigreeReps;

use libsvm::{SvmInit, SvmTrainer, KernelInit, ModelInit};

mod error;
use error::{SVMError, SvmBuilderError};

pub struct LibSvmBuilder<'a> {
    inner      : SvmTrainer,
    features   : Option<Vec<[f64; 1]>>,
    labels     : Option<Vec<f64>>,
    label_order: Option<&'a [&'a String]>,
    mu         : f64,
    sigma      : f64,
}

impl<'a> Default for LibSvmBuilder<'a> {
    fn default() -> Self {
        // ---- Initialize binary SVM with linear kernel.
        let svm_initializer  = SvmInit {
            kernel                : Some(KernelInit::Linear),
            model                 : Some(ModelInit::CSvc { cost: Some(10.) }),
            label_weights         : None,
            termination_eps       : Some(0.001),
            probability_estimates : Some(true),
            shrinking             : Some(true),
            cache_size            : Some(2048)
        };
        let svm_trainer = svm_initializer.build()
            .expect("Failed to initialize libsvm SVMTrainer");
        Self { inner: svm_trainer, features: None, labels: None, label_order: None, mu: 0.0, sigma: 1.0}
    }
}

impl<'a> LibSvmBuilder<'a> {
    pub fn sims(&mut self, pedigrees: &PedigreeReps) -> &mut Self {
        let mut data: Vec<[f64; 1]> = Vec::with_capacity(pedigrees.len());
        data.extend(
            pedigrees.iter()
            .flat_map(|ped| ped.comparisons.iter()
                .map(|cmp| [cmp.get_avg_pwd(); 1])
            )
        );
        self.features = Some(data);
        self
    }

    pub fn label_order(&mut self, label_order: &'a [&String]) -> &mut Self {
        self.label_order = Some(label_order);
        self
    }

    pub fn labels(&mut self, pedigrees: &PedigreeReps, label_treshold: &usize) -> &mut Self {
        let Some(label_order) = self.label_order else {
            panic!("No underflying label order found.")
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

    /// Standardize vector (mu=0.0, sigma=1.0)
    pub fn scale(&mut self) -> Result<&mut Self> {
        let Some(ref mut features) = self.features else {
            return Err(SvmBuilderError::MissingFeatures).loc("While scaling SVM features")
        };
        let len   = features.len() as f64;
        let mu    = features.iter().map(|v| v[0]).sum::<f64>() / len ;
        let sigma = (features.iter().map(|v| (v[0] - mu).powf(2.0)).sum::<f64>() / len).sqrt();

        trace!("Scaled SVM features with:  mu={mu} and sigma={sigma}");
        if sigma > 0.0 {
            features.iter_mut().for_each(|v| v[0] = (v[0] - mu)/sigma);
            (self.mu, self.sigma) = (mu, sigma);
        } else {
            warn!("Invalid standard deviation while applying SVM feature scaling: {sigma}")
        }

        Ok(self)
    }

    pub fn build(&self) -> Result<LibSvm> {
        let loc_msg = "While attempting to build SVM from builder.";
        let Some(ref labels) = self.labels else {
            return Err(SvmBuilderError::MissingLabels).loc(loc_msg)
        };

        let Some(ref features) = self.features else {
            return Err(SvmBuilderError::MissingFeatures).loc(loc_msg)
        };

        trace!("Fitting SVMPredictor...");
        let svm = self.inner.fit(&features[..], &labels[..])
            .map_err(SvmBuilderError::FitSvm)
            .loc(loc_msg)?;
        
        Ok(LibSvm { inner: svm, mu: self.mu, sigma: self.sigma })
    }
}
pub struct LibSvm {
    inner: libsvm::SvmPredictor,
    mu   : f64,
    sigma: f64,
}

impl LibSvm {
    /// Directly create an SVM from a pedigree vector.
    #[allow(unused)]
    pub fn from_comparisons(pedigrees: &PedigreeReps, label_treshold: &usize, label_order: &[&String]) -> Result<Self> {

        let svm = LibSvmBuilder::default()
            .sims(pedigrees)
            .label_order(label_order)
            .labels(pedigrees, label_treshold)
            .build()
            .loc("While attempting to build SVM.")?;

        Ok( Self{inner: svm.inner, mu: svm.mu, sigma: svm.sigma} )
    }


    pub fn predict_probs(&mut self, avg: f64) -> Result<Vec<(f64, Vec<f64>)>> {
        let x: [f64; 1] = [(avg - self.mu) / self.sigma];
        let x: &[f64] = x.as_ref();

        self.inner.predict_with_probability(x)
            .map_err(SVMError::PredictProbs)
            .with_loc(|| format!("While attempting to predict value {avg}"))
    }

    #[allow(unused)]
    pub fn predict(&self, avg: f64) -> Result<f64>{
        let x: [f64; 1 ] = [avg];
        let x: &[f64] = x.as_ref();
        Ok(self.inner.predict(x)
            .map_err(SVMError::PredictBool)
            .with_loc(|| format!("While attempting to predict value {avg}"))?[0]
        )   
    }

    pub fn true_label_idx(&self) -> Result<usize, SVMError> {
        self.inner.labels().last().map(|idx| *idx as usize).ok_or(SVMError::InvalidLabels)
    }
}

