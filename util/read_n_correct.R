# Custom function for upset plots
#
# 2020
# CCDL - C. Savonen

## Set up special functions

read_n_correct <- function(filename,
                           wavenum_min,
                           wavenum_max,
                           peak_threshold = 3,
                           min_cutoff = NA, 
                           max_cutoff = NA,
                           max_scale = NA, 
                           ridge_length = 32,
                           window_size = 500,
                           return_plots = FALSE, 
                           more_peak_info = TRUE) {
  # Will background correct and identify peaks for a given vibration dataset
  # saved as a CSV file from
  #
  # Args:
  #   filename : A file path to to a vibration CSV file that has two columns
  #             The first column is wavenumber and the second is intensity.
  #   wavenum_min/max : two numbers which specify the range of numbers to keep in
  #                    the dataset. Wavenumbers outside these ranges will be
  #                    removed from the data.frame
  #   min/max_cutoff: If you'd like the data to be filtered before correction then 
  #                   provide the ranges of wavelengths you'd like to be included.
  #   ridge_length: passed directly to identifyMajorPeaks() function
  #   return_plots: TRUE/FALSE should the plots be returned?
  #
  # Returns:
  # A data.frame with background corrected intensities and peaks identified.
  # Get file.path
  #
  
  filename <- file.path(filename)

  # Check if the file actually exists
  if(!file.exists(filename)){
    stop(paste("No file at the given file path."))
  }

  # Read in the CSV data file colnames will be ignored
  vibration_df <- readr::read_csv(filename, col_names = FALSE) %>%
    # Save the sensible column names
    dplyr::rename("wavenumber" = X1,
                  "intensity" = X2)
  
  # Filter data based on optional min_cutoff
  if (!is.na(min_cutoff)){
  vibration_df <- vibration_df %>% 
    dplyr::filter(wavenumber > min_cutoff)
  }
  
  # Filter data based on optional max_cutoff
  if (!is.na(max_cutoff)){
    vibration_df <- vibration_df %>% 
      dplyr::filter(wavenumber < max_cutoff)
  }
  
  # If no max_scale is specified, try it out 
  if (is.na(max_scale)) {
    # Establish the scale
    scales <- seq(1, 120, 1)
  
    # Try cwt with 120 and see if it fails
    try_it <- try(cwt(vibration_df$intensity,
                      # The scales at which to perform CWT
                      scales = scales,
                      # Which wavelet base to use
                      wavelet = 'mexh'), 
                  silent = TRUE)
  
    # Declare max_scale based on the width
    max_scale <- ifelse(grepl("Error", try_it)[1], 63, 120)
  }
  
  # Print out scale used
  print(paste("max_scale set to:", max_scale))
  
  # Establish the scale
  scales <- seq(1, max_scale, 1)

  # Find Continuous Wavelet Transform
  w_coefs <-  cwt(vibration_df$intensity,
                  # The scales at which to perform CWT
                  scales = scales,
                  # Which wavelet base to use
                  wavelet = 'mexh')

  # Identify the local maximum of each column in 2-D CWT coefficients matrix
	# by using a slide window
	local_max <- getLocalMaximumCWT(w_coefs,
	                                # Minimum window
	                                minWinSize = 0.5)

  # Identify ridges by connecting the local maximum of 2-D
	ridge_list <- getRidge(local_max,
	                       gapTh = 3,
	                       skip = NULL)
  if (return_plots) {
    # Make this dat into an image
	   print(
	  image(1:nrow(w_coefs),
	        scales,
	        w_coefs,
	        col = rainbow(256),
	        axes = FALSE,
	        xlab = 'index',
	        ylab = 'CWT coefficient scale',
	        main = 'CWT coefficients')
	  )

	  # Plot the ridge list
	  plotRidgeList(ridge_list)
  }
  # Identify the major peaks
	major_peak_info <- identifyMajorPeaks(vibration_df$intensity,
	                                      ridge_list,
	                                      w_coefs,
	                                      SNR.Th = peak_threshold,
	                                      ridgeLength = ridge_length,
	                                      nearbyWinSize = window_size,
	                                      nearbyPeak = TRUE)
  # Estimate peak width
	peak_width <- widthEstimationCWT(vibration_df$intensity,
	                                 major_peak_info)

	# Make the data.frame with the corrected data columns.
	vibration_df <- vibration_df %>%
	  dplyr::mutate(
	    # Make background a new column based on intensity column
	    background = baselineCorrectionCWT(
	      intensity,
	      peak_width,
	      lambda = 1000,
	      differences = 1),
	    # Save corrected intensity as a new column
	    corrected_intensity = intensity - background,
	    # Save binned column
	    binned = lowess(wavenumber,
	                    corrected_intensity,
	                    f = 0.03)$x,
	    # Save uncorrected bin column
	    uncorrected_bin = lowess(wavenumber,
	                             intensity,
	                             f = 0.03)$x,
	    # Save new column that says which are major peaks
	    major_peak = (1:nrow(.) %in% major_peak_info$peakIndex),
	    peak_index = 1:nrow(.),
	    wavenumber = as.numeric(wavenumber)
	    ) %>%
	  dplyr::filter(wavenumber > wavenum_min,
	                wavenumber < wavenum_max) 
	  
	  # If more_peak_info is TRUE
	  if (more_peak_info) {
	    
	    # Extract SNR based on peak index
	    peak_snr_df <- data.frame(
	      peak_index = major_peak_info$peakIndex, 
	      peak_snr = major_peak_info$peakSNR[match(names(major_peak_info$peakIndex), names(major_peak_info$peakSNR))]
	      )
	    
	    # Join it on
	    vibration_df <- vibration_df %>% 
	     dplyr::left_join(peak_snr_df, by = "peak_index")
	  }
	
	
	# Return the corrected and filtered dataset
	return(vibration_df)
}
