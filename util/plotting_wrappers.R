# Custom function for upset plots
#
# 2020
# CCDL - C. Savonen

plot_ind_peaks <- function(vibration_df) {
  # Given a vibration data.frame that has been corrected already using
  # `read_n_correct` function, make an individual peak shape for ultimate plot
  #
  # Args:
  #   vibration_df : A corrected vibration dataset as a data.frame with
  #                  `corrected_intensity`, `wavenumber` and `major_peak` as
  #                  columns.
  #   x_lab: x axis label
  #   y_lab: y axis label
  #
  # Returns:
  # A data.frame with background corrected intensities.

  # Make a plot with no axises or other stuff, but has the peaks
	ggplot2::ggplot(vibration_df) +
	  ggplot2::geom_line(ggplot2::aes(x = wavenumber,
	                                  y = corrected_intensity)) +
	  ggplot2::theme_classic() +
    ggplot2::xlim(c(round(min(vibration_df$wavenumber)),
                    round(max(vibration_df$wavenumber)))) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    # Make everything blank because we just want the peak shape
    ggplot2::theme(legend.position = "none",
                   # Make plot margins smaller
                   plot.margin = grid::unit(c(0, 0, 0, 0), "cm"),
                   axis.line.y = ggplot2::element_blank(),
                   # set transparency
                   panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   plot.background =ggplot2::element_rect(fill = "transparent", colour = NA)
                   )
}

remove_x_axis <- function(gplot) {
  gplot <- gplot + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  axis.line.x = ggplot2::element_blank())

  return(gplot)
}

remove_y_axis <- function(gplot) {
  gplot <- gplot + ggplot2::theme(
  axis.title.y = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank())

  return(gplot)
}

remove_axes <- function(gplot) {
  gplot <- gplot + ggplot2::theme(
  axis.title.y = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank(),
  axis.title.x = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  axis.line.x = ggplot2::element_blank())

  return(gplot)
}

add_text_labels <- function(gplot, text) {
  gplot <- gplot +
    ggplot2::annotate("text",
                      label = text,
                      x = max(gplot$data$wavenumber) - 150,
                      y = max(gplot$data$corrected_intensity) - 100,
                      size = 5,
                      colour = "red",
                      fontface = "bold")
  return(gplot)
}

assemble_main_plot <- function(peak_plots_list, 
                               row_pressures,
                               main_title = "The background-correction of Raman Spectra",
                               x_lab = "Wavelength / nm",
                               y_lab = "Raman Intensity/Arbitr. Units") {
  # Given a list of plots (output from `plot_ind_peaks`) and a vector of pressures
  # that are in the names of that plot, assemble the main plot 
  #
  # Args:
  #   peak_plots: a list of individual plots with the pressures in the names
  #   row_pressure: a vector of pressures that will be used for labeling and for 
  #                 assembling each row. Each row will be one of the pressure sets.
  #   main_title, x_lab, y_lab: character strings to be used as those labels respectively. 
  #                             There are defaults for each. 
  #
  # Returns:
  # A plot with all the individual peak plots assembled into one
  
  # Take all the row pressures and sort them in decreasing order. 
  row_pressures <- as.numeric(unique(row_pressures))
  row_pressures <- sort(row_pressures, decreasing = TRUE)
  
  # Establish which pressure is the smallest
  min_pressure <- min(row_pressures)
  
  # Assemble each row based on its pressure noted in the name
  each_pressure_plots <- 
    lapply(row_pressures, function(a_row){
      
      # Use the names in the list to find which ones have the pressures 
      row_plot_names <- grep(a_row, names(peak_plots_list), value = TRUE)
      
      row_plots <- peak_plots_list[row_plot_names]
      
      # If its the smallest pressure, we still want the x axes
      if (a_row == min_pressure){
        # Get the first plot set up 
        the_row <- remove_y_axis(row_plots[[1]])
        
        # Add on the subsequent plots
        for (plot in 2:length(row_plots)){
          the_row <- the_row + remove_y_axis(row_plots[[plot]])
        }
      # All other pressure rows we want to remove both y and x axes
      } else {
        # Get the first plot set up
        the_row <- remove_axes(row_plots[[1]])
        
        # Add on the subsequent plots
        for (plot in 2:length(row_plots)){
          the_row <- the_row + remove_axes(row_plots[[plot]])
        }
      }
      
      # Add the pressure label here
      the_row <- add_text_labels(the_row, text = paste(a_row, "GPa"))
      
      return(the_row)
    })
  
  # Now assemble the main plot with each row's plot
  ultimate_plot <- each_pressure_plots[[1]]
  
  for (row in 2:length(each_pressure_plots)){
    ultimate_plot <- ultimate_plot / each_pressure_plots[[row]]
  }
  
  # Set up x label
  x_lab <- wrap_elements(grid::textGrob(x_lab)) + 
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(.001, .001, .001, .001), "cm")
    )
  
  # Set up y label
  y_lab <- wrap_elements(grid::textGrob(y_lab, rot = 90)) + 
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(.01, .01, .01, .01), "cm")
    )
  
  # Tack on the labels to the main plot
  ultimate_plot <- 
    # Start with y label
    (y_lab + ultimate_plot) +
    plot_layout(
      widths = c(0.1, 5) # Adjust spacing 
    ) 
  
  ultimate_plot <- 
    (ultimate_plot / x_lab) + 
    # Adjust spacing 
    plot_layout(
      heights = c(5, 0.1) # Adjust spacing 
    )
  
  ultimate_plot <- ultimate_plot + 
    plot_annotation(
      title = ggplot2::element_text(main_title),
    )
  
  return(ultimate_plot)
}
