use ndarray::{Array1};
use crate::pedigrees::PedigreeReps;

use smartcore::{svm::{LinearKernel, svc::{SVC, SVCParameters}}, linalg::basic::matrix::DenseMatrix};

/// Benchmark with smartcore library.
/// Extremely slow. Am I doing something wrong ???
#[allow(unused)]
pub fn smartcore_svm(pedigrees: &PedigreeReps, label_treshold: &usize, label_order: &[&String], preds: &[f64]) {
    let data: Vec<Vec<f64>> = pedigrees.iter()
        .flat_map(|ped| ped.comparisons.iter()
            .map(|cmp| vec![cmp.get_avg_pwd()])
        )
        .collect::<Vec<_>>();

    let data = DenseMatrix::from_2d_vec(&data);
    
    let targets: Array1<i8> = Array1::from_iter(pedigrees.iter()
        .flat_map(|ped| ped.comparisons.iter()
            .map(|cmp| {
                let label_idx = label_order.iter()
                    .position(|lab| *lab == &cmp.label )
                    .expect("Invalid label index while building smartcore SVM");
                match label_idx > *label_treshold {
                    true   => 1,
                    false => -1
                }
            })
        )
    );
    
    let parameters = SVCParameters::default()
        .with_kernel(LinearKernel)
        .with_c(1.0);

    println!("Fit...");
    let svm = SVC::fit(&data, &targets, &parameters)
        .expect("Failed to fit Smartcore SVC");
    let pred  = svm.predict(&DenseMatrix::from_2d_array(&[preds]));
    let probs = svm.decision_function(&DenseMatrix::from_2d_array(&[preds]))
        .expect("Failed to extract decision values");
    println!("Done. Decisions: {probs:?}");
}
