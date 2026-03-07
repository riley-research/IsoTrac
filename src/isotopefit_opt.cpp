#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double isotopefit_opt_cpp(const DataFrame& spectrum,
                          const DataFrame& isotope_template,
                          double poi,
                          double ppm_tolerance,
                          double cosine_score_correction,
                          NumericVector minimal_isotope_sequence,
                          std::string score) {

  // 1. Extract and Clone (Safe modification)
  NumericVector spec_mz = spectrum["mz"];
  NumericVector spec_int = spectrum["centroided_intensity"];
  NumericVector temp_mz = clone(as<NumericVector>(isotope_template["mz"]));
  NumericVector temp_int = clone(as<NumericVector>(isotope_template["intensity"]));
  NumericVector temp_seq = isotope_template["isotope_seq"];
  NumericVector temp_pct = isotope_template["percent"];

  int n_temp = temp_mz.size();
  int n_spec = spec_mz.size();

  // --- Step 1: Find Anchors & Calculate Shifts ---
  double max_t_int = -1.0;
  double mz_at_max = 0.0;
  double intensity_at_100pct = 0.0;

  for(int i = 0; i < n_temp; ++i) {
    if(temp_int[i] > max_t_int) {
      max_t_int = temp_int[i];
      mz_at_max = temp_mz[i];
    }
    if(std::abs(temp_pct[i] - 100.0) < 1e-6) {
      intensity_at_100pct = temp_int[i];
    }
  }
  if(intensity_at_100pct <= 0) intensity_at_100pct = max_t_int;

  // Calculate m/z shift
  double difmz = poi - mz_at_max;

  // Robust lookup for max_spec_int at the POI
  double max_spec_int = 0.0;
  auto it_poi = std::lower_bound(spec_mz.begin(), spec_mz.end(), poi);

  // Check the found index and the one before it for the closest match within tolerance
  int poi_idx = std::distance(spec_mz.begin(), it_poi);
  double best_dist = 1e9;

  for(int j = poi_idx - 1; j <= poi_idx; ++j) {
    if(j >= 0 && j < n_spec) {
      double d = std::abs(spec_mz[j] - poi);
      if(d < best_dist && (d / poi * 1e6) <= ppm_tolerance) {
        best_dist = d;
        max_spec_int = spec_int[j];
      }
    }
  }

  // Calculate Intensity Scaling
  double difint = (intensity_at_100pct > 0) ? (max_spec_int / intensity_at_100pct) : 1.0;

  // --- NEW: Global Update (Applied ONCE here) ---
  for(int i = 0; i < n_temp; ++i) {
    temp_mz[i] += difmz;
    temp_int[i] *= difint;
  }

  // --- Step 2 & 3: Match Peaks using UPDATED temp_mz ---
  std::vector<double> matched_temp_int;
  std::vector<double> matched_spec_int;
  std::vector<int> found_seq;

  for(int i = 0; i < n_temp; ++i) {
    double target_mz = temp_mz[i]; // No + difmz here, it's already in temp_mz

    auto it = std::lower_bound(spec_mz.begin(), spec_mz.end(), target_mz);
    int idx = std::distance(spec_mz.begin(), it);

    int closest = -1;
    double min_d = 1e9;

    // Check idx and idx-1 for closest m/z
    for(int j = idx - 1; j <= idx; ++j) {
      if(j >= 0 && j < n_spec) {
        double d = std::abs(spec_mz[j] - target_mz);
        if(d < min_d) {
          min_d = d;
          closest = j;
        }
      }
    }

    if(closest != -1) {
      double ppm_err = (min_d / target_mz) * 1e6;
      if(ppm_err <= ppm_tolerance) {
        matched_temp_int.push_back(temp_int[i]);
        matched_spec_int.push_back(spec_int[closest]);
        found_seq.push_back((int)temp_seq[i]);
      }
    }
  }

  // --- Step 4: Validate Minimal Sequence ---
  for(int i = 0; i < minimal_isotope_sequence.size(); ++i) {
    if(std::find(found_seq.begin(), found_seq.end(), (int)minimal_isotope_sequence[i]) == found_seq.end()) {
      return 0.0;
    }
  }

  // Identify monoisotopic (0) and expand (Contiguous check)
  std::vector<double> final_temp, final_spec;
  int zero_idx_in_template = -1;
  for(int i = 0; i < n_temp; ++i) { if(temp_seq[i] == 0) { zero_idx_in_template = i; break; } }
  if(zero_idx_in_template == -1) zero_idx_in_template = 0;

  // Forward
  for(int i = zero_idx_in_template; i < n_temp; ++i) {
    auto it = std::find(found_seq.begin(), found_seq.end(), (int)temp_seq[i]);
    if(it != found_seq.end()) {
      int d = std::distance(found_seq.begin(), it);
      final_temp.push_back(matched_temp_int[d]);
      final_spec.push_back(matched_spec_int[d]);
    } else break;
  }
  // Backward
  for(int i = zero_idx_in_template - 1; i >= 0; --i) {
    auto it = std::find(found_seq.begin(), found_seq.end(), (int)temp_seq[i]);
    if(it != found_seq.end()) {
      int d = std::distance(found_seq.begin(), it);
      final_temp.insert(final_temp.begin(), matched_temp_int[d]);
      final_spec.insert(final_spec.begin(), matched_spec_int[d]);
    } else break;
  }

  if(final_temp.size() < 5) return 0.0;

  // --- Step 5: Scoring (Pearson / Cosine) ---
  if(score == "pearson") {
    double sum_a = 0, sum_b = 0;
    for(size_t i=0; i<final_temp.size(); ++i) { sum_a += final_temp[i]; sum_b += final_spec[i]; }
    double m_a = sum_a / final_temp.size(), m_b = sum_b / final_temp.size();
    double num = 0, den_a = 0, den_b = 0;
    for(size_t i=0; i<final_temp.size(); ++i) {
      num += (final_temp[i]-m_a)*(final_spec[i]-m_b);
      den_a += std::pow(final_temp[i]-m_a, 2);
      den_b += std::pow(final_spec[i]-m_b, 2);
    }
    return (den_a > 1e-12 && den_b > 1e-12) ? num / std::sqrt(den_a * den_b) : 0.0;

  } else if(score == "cosine") {
    double min_val = final_temp[0];
    for(size_t i=0; i<final_temp.size(); ++i) {
      if(final_temp[i] < min_val) min_val = final_temp[i];
      if(final_spec[i] < min_val) min_val = final_spec[i];
    }
    double offset = min_val * cosine_score_correction;
    double dot = 0, mag_a = 0, mag_b = 0;
    for(size_t i=0; i<final_temp.size(); ++i) {
      double a = final_temp[i] - offset;
      double b = final_spec[i] - offset;
      dot += a * b; mag_a += a * a; mag_b += b * b;
    }
    return (mag_a > 1e-12 && mag_b > 1e-12) ? dot / (std::sqrt(mag_a) * std::sqrt(mag_b)) : 0.0;
  }

  return 0.0;
}
