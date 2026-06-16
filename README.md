# Installation 

Run the following code to install `IsoTrac`:

 ```R 
# Install pak if you don't have it yet if
 (!requireNamespace("pak", quietly = TRUE)) { install.packages("pak") } 

# Install IsoTrac from the GitHub repository 
pak::pkg_install("riley-research/IsoTrac")
```

# Input files
Preparing the input file is most easily done by using [PTsliCR](https://github.com/riley-research/PTsliCR). PTsliCR outputs a file `PTCR_cleaned_for_extraction.txt` that is directly compatible with IsoTrac for deconvolution of the masses.

# Usage

An example of how to use IsoTrac:

 ```R 
library(IsoTrac)

aa_seq <- toupper("qeecvcenyklavncfvnnnrqcqctsvgaqntvicsklaakclvmkaemngsklgrrakpegalqnndglydpdcdesglfkakqcngtstcwcvntagvrrtdkdteitcservrtywiiielkhkarekpydskslrtalqkeittryqldpkfitsilyennvitidlvqnssqktqndvdiadvayyfekdvkgeslfhskkmdltvngeqldldpgqtliyyvdekapefsmqglkhhhhhh")

protein_mass <- 28211.86

mgfPath <- "C:/PTCR_cleaned_for_extraction.txt"

my_mgf <- import_spectra(mgfPath, intensity_threshold = 50, charge_range = 1:16,
               scan_range = c(1,20))

plot_spectra(my_mgf, "C:/Analysis/")

my_fdr <- generate_fdr_spectra(my_mgf, isotope_peaks_included = 4)

threshold <- get_fdr_threshold(my_mgf, my_fdr, aa_seq, protein_mass,
                               minimal_isotope_sequence = c(-3, -2, -1, 0, 1, 2, 3), isotope_peaks_included = 4,
                               number_of_cores = 2, score = "pearson", plot_bins = 0.06)


results <- extract_isotopes(my_mgf, aa_seq, protein_mass, threshold = 0.9, score = "pearson", number_of_cores =  2,
                                 minimal_isotope_sequence = c(-3, -2, -1, 0, 1, 2, 3), isotope_peaks_included = 11, nrmse_cutoff = 45,
                                 output_folder = "C:/Analysis")

results$founddf
results$nonefounddf
```

For further documentation, please see the vignettes and manual in the XXX folder.

The documentation is also available through R by running ? followed by the function name. For example:

 ```R 
?extract_isotopes
```

# Output

The following will be output to the user-supplied output folder of the extract_isotopes function:

- `matched_peaks.csv`: all data, including correlation scores and intact masses, of the matched isotope distributions.
- `matched_spectra`: the theoretical and imperical isotope cluster overlaps of all peaks in the matched_peaks.csv.
- `unmatched_peaks`: the data for slices without any matches.
- `unmatched_spectra`: the theoretical and imperical isotope cluster overlaps of all peaks in the unmatched_peaks.csv.
- `separated_matched_spectra` and `separated_unmatched_spectra` folders: the same visualizations as the matched_peaks and unmatched_peaks, but in separate pdf files.

The extract_isotopes function will also return a list in R that holds the `founddf` and `nonefounddf` dataframes. These correspond to the matched_peaks.csv and unmatched_peaks.csv file.

# Questions, feedback, missing capabilities

Please use the Issues tab on GitHub or send an email to tveth@uw.edu
