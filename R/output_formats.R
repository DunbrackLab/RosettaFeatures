# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#' @export
get_output_formats <- function(requested_output_formats, add_footer){
	output_formats %>%
		dplyr::semi_join(dplyr::data_frame(id=requested_output_formats), by="id") %>%
		dplyr::mutate(
			add_footer = add_footer & accepts_footer)
}

