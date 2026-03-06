#include <Rcpp.h>
#include <vector>
#include <string> // Add this

using namespace Rcpp;

// [[Rcpp::export]]
double isotopefit_opt_cpp(const DataFrame& spectrum,
                      const DataFrame& isotope_template,
                      double poi,
                      double ppm_tolerance,
                      double cosine_score_correction,
                      NumericVector minimal_isotope_sequence,
                      std::string score) { // Use std::string here

  // Extract columns (use as<> to be explicit and safe)
  NumericVector spec_mz = spectrum["mz"];
  NumericVector spec_int = spectrum["centroided_intensity"];

  // Use clone() only for vectors you intend to change (like the template mz)
  NumericVector temp_mz = clone(as<NumericVector>(isotope_template["mz"]));
  NumericVector temp_int = clone(as<NumericVector>(isotope_template["intensity"]));
  NumericVector temp_seq = isotope_template["isotope_seq"];
  NumericVector temp_pct = isotope_template["percent"];

  int n_temp = temp_mz.size();
  int n_spec = spec_mz.size();

  // --- Step 1: Align template to POI & Calculate Scaling ---
  double max_t_int = -1.0;
  double mz_at_max = 0.0;
  double intensity_at_100pct = 0.0;

  // First pass: Find anchors
  for(int i = 0; i < n_temp; ++i) {
    if(temp_int[i] > max_t_int) {
      max_t_int = temp_int[i];
      mz_at_max = temp_mz[i];
    }
    // Match R: find intensity where percent is 100
    if(std::abs(temp_pct[i] - 100.0) < 1e-6) {
      intensity_at_100pct = temp_int[i];
    }
  }

  // If intensity_at_100pct wasn't found by percent, fallback to max theoretical intensity
  if(intensity_at_100pct <= 0) intensity_at_100pct = max_t_int;

  // Calculate shifts (Matches R exactly)
  double difmz = poi - mz_at_max;

  double max_spec_int = 0.0;
  for(int i = 0; i < n_spec; ++i) {
    if(std::abs(spec_mz[i] - poi) < 1e-6) {
      max_spec_int = spec_int[i];
      break;
    }
  }

  double difint = (intensity_at_100pct > 0) ? (max_spec_int / intensity_at_100pct) : 1.0;

  // --- NEW: Global Update (Matches R: isotope_template$mz <- isotope_template$mz + difmz) ---
  for(int i = 0; i < n_temp; ++i) {
    temp_mz[i] += difmz;
    temp_int[i] *= difint;
  }

  // --- Step 2 & 3: Match Peaks and Filter by PPM ---
  std::vector<double> matched_temp_int;
  std::vector<double> matched_spec_int;
  std::vector<int> found_seq;

  for(int i = 0; i < n_temp; ++i) {
    double shifted_mz = temp_mz[i] + difmz;

    // Binary search for closest m/z
    auto it = std::lower_bound(spec_mz.begin(), spec_mz.end(), shifted_mz);
    int idx = std::distance(spec_mz.begin(), it);

    int closest = -1;
    if (idx < n_spec) closest = idx;
    if (idx > 0 && (closest == -1 || std::abs(spec_mz[idx-1] - shifted_mz) < std::abs(spec_mz[idx] - shifted_mz))) {
      closest = idx - 1;
    }

    if(closest != -1) {
      double s_mz = spec_mz[closest];
      double ppm_err = (shifted_mz - s_mz) / s_mz * 1e6;
      if(std::abs(ppm_err) <= ppm_tolerance) {
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

  // Trim missing sequence (Contiguous check)
  // We identify the monoisotopic (0) and expand outwards
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

  // --- Step 5: Scoring ---
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
    return (den_a > 0 && den_b > 0) ? num / std::sqrt(den_a * den_b) : 0.0;

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
    return (mag_a > 0 && mag_b > 0) ? dot / (std::sqrt(mag_a) * std::sqrt(mag_b)) : 0.0;
  }

  return 0.0;
}
