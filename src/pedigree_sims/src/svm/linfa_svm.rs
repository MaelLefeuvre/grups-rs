#![allow(unused)]
use linfa::{prelude::*, Platt, platt_scaling::{platt_predict, PlattValidParams}, MultiClassModel};
use linfa_svm::{Svm, SvmParams};

use ndarray::{Array1, Array2};
use crate::pedigrees::PedigreeReps;


/// Benchmark with linfa library.
/// Slowish + can't figure a way to obtain decent Platt scaling.
/// That's a shame, cause this library seems more robust, but the documentation is quite lacking...
#[allow(unused)]
pub struct LinfaSVM {
    inner: Svm<f64, bool>
    //inner: Platt<f64, Svm<f64, f64>>
}

#[allow(unused)]
impl LinfaSVM {
    pub fn from_comparisons(pedigrees: &PedigreeReps, label_treshold: &usize, label_order: &[&String]) -> Self {
        
        let data = pedigrees.iter()
            .flat_map(|ped| ped.comparisons.iter()
                .map(|cmp| [cmp.get_avg_pwd()])
            )
            .collect::<Vec<_>>()
            .into();

        let targets = Array1::from_iter(
            pedigrees.iter()
                .flat_map(|ped| ped.comparisons.iter()
                    .map(|cmp| 
                        label_order.iter()
                            .position(|lab| *lab == &cmp.label )
                            .expect("Invalid Label Found while building SVM labels.")
                    )
                )
                .collect::<Vec<_>>()
            );

        Self::new(data, targets, label_treshold)
    }

    pub fn new(data: Array2<f64>, targets: Array1<usize>, label_treshold: &usize) -> Self {
        let feature_names = ["avg", "label"];
        let (train, _) = Dataset::new(data, targets)
            //.with_feature_names(feature_names)
            .map_targets(|label| label > label_treshold)
            .split_with_ratio(1.0);

        println!(
            "Fit SVM classifier with #{} training points",
            train.nsamples()
        );

        // fit a SVM with C value 7 and 0.6 for positive and negative classes
        let svm = Svm::<f64, bool>::params()
            //.shrinking(true) // Faster, but less accurate
            //.pos_neg_weights(50000., 5000.)
            //.with_platt_params(Platt::params())
            .linear_kernel()
            .fit(&train)
            .expect("Failed to fit linfa SVM.");     

        //Platt::predict(svm, ndarray::array![0.13]);
        //let model = Platt::params()
        //    .fit_with(svm.clone(), &train.map_targets(|t| *t > 0.5))
        //    .expect("Failed to fit linfa SVM with Platt scaling");;
        Self{inner: svm}
    }

    pub fn predict(&self, avg: f64) -> bool {
        let arr = ndarray::array![avg];
        self.inner.predict(arr)
    }

}
