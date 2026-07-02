use scicdh::statistics::DataSet;

#[test]
fn test_mean_calculation() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let dataset = DataSet::new(data);
    
    let mean = dataset.mean().expect("Failed to calculate mean");
    assert_eq!(mean, 3.0);
}

#[test]
fn test_median_odd_length() {
    let data = vec![7.0, 1.0, 3.0, 9.0, 5.0]; // Unsorted input
    let dataset = DataSet::new(data);
    
    let median = dataset.median().expect("Failed to calculate median");
    assert_eq!(median, 5.0); // Sorted middle: 1, 3, [5], 7, 9
}

#[test]
fn test_median_even_length() {
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    let dataset = DataSet::new(data);
    
    let median = dataset.median().expect("Failed to calculate median");
    assert_eq!(median, 4.5); // (4.0 + 5.0) / 2.0
}

#[test]
fn test_mode_unimodal() {
    let data = vec![1.0, 2.0, 2.0, 3.0, 4.0];
    let dataset = DataSet::new(data);
    
    let mode = dataset.mode().expect("Failed to calculate mode");
    assert_eq!(mode, Some(2.0));
}

#[test]
fn test_mode_uniform_distribution() {
    // A uniform sequence has no distinct majority frequency element
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let dataset = DataSet::new(data);
    
    let mode = dataset.mode().expect("Failed to calculate mode");
    assert_eq!(mode, None);
}

#[test]
fn test_range_calculation() {
    let data = vec![-5.0, 2.0, 10.0, 0.5];
    let dataset = DataSet::new(data);
    
    let range = dataset.range().expect("Failed to calculate range");
    assert_eq!(range, 15.0); // 10.0 - (-5.0)
}

#[test]
fn test_sample_variance_and_std_deviation() {
    let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    let dataset = DataSet::new(data);
    
    let variance = dataset.sample_variance().expect("Failed variance evaluation");
    let std_dev = dataset.sample_std_deviation().expect("Failed std dev evaluation");
    
    // Unbiased sample variance for this set is exactly 4.571428571428571
    assert!((variance - 4.571428).abs() < 1e-5);
    assert!((std_dev - 4.571428_f64.sqrt()).abs() < 1e-5);
}

#[test]
fn test_skewness_and_excess_kurtosis_flat_uniform() {
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    let dataset = DataSet::new(data);
    
    let skewness = dataset.skewness().expect("Failed skewness");
    let kurtosis = dataset.kurtosis().expect("Failed kurtosis");
    
    // Symmetric data yields exactly 0 skewness
    assert_eq!(skewness, 0.0);
    // Flat uniform distribution yields exactly -1.2242424242424244 excess kurtosis
    assert!((kurtosis - (-1.224242)).abs() < 1e-5);
}

#[test]
fn test_pearson_correlation() {
    let dataset_x = DataSet::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    let y_vector = vec![2.0, 4.0, 6.0, 8.0, 10.0]; // Perfectly linearly correlated (y = 2x)
    
    let r = dataset_x.pearson_correlation_coeff(&y_vector).expect("Failed correlation");
    assert!((r - 1.0).abs() < 1e-7);
}

#[test]
fn test_has_and_into_set() {
    let data = vec![1.0, 1.0, 2.0, 3.0, 3.0, 3.0];
    let dataset = DataSet::new(data);
    
    assert!(dataset.has(2.0, None).unwrap());
    assert!(!dataset.has(5.0, None).unwrap());
    
    let unique_set = dataset.into_set().expect("Failed deduplication");
    assert_eq!(unique_set.len(), 3);
    assert!(unique_set.contains(&1.0));
    assert!(unique_set.contains(&2.0));
    assert!(unique_set.contains(&3.0));
}

#[test]
fn test_empty_dataset_error_handling() {
    let dataset = DataSet::new(vec![]);
    
    assert!(dataset.mean().is_err());
    assert!(dataset.median().is_err());
    assert!(dataset.mode().is_err());
    assert!(dataset.range().is_err());
    assert!(dataset.sample_variance().is_err());
}

