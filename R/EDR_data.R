#' Ecological Dynamic Regime data
#'
#' @description
#' Example datasets to characterize and compare EDRs, including abundance data,
#' state, segment, and trajectory dissimilarity matrices for 90 artificial communities
#' belonging to three different EDRs.
#'
#' @format
#' List of three nested sublists (`"EDR1"`, `"EDR2"`, and `"EDR3"`), each
#' associated with one EDR, including the following elements:
#' * `abundance`: Data table with 15 columns and one row for each community state:
#'     + `EDR`: Integer indicating the identifier of the EDR.
#'     + `traj`: Integer containing the identifier of the trajectory for each
#'     artificial community in the corresponding EDR. Each trajectory represents
#'     a different sampling unit.
#'     + `state`: Integer indicating the observations or states of each community.
#'     The sequence of states of a given community forms a trajectory.
#'     + `sp1, ..., sp12`: Vectors containing species abundances for each community
#'     state.
#' * `state_dissim`: Object of class `dist` containing Bray-Curtis dissimilarities
#' between every pair of states in `abundance` (see Details).
#' * `segment_dissim`: Object of class `dist` containing the dissimilarities
#' between every pair of trajectory segments in `abundance` (see Details).
#' * `traj_dissim`: Object of class `dist` containing the dissimilarities
#' between every pair of community trajectories in `abundance` (see Details).
#'
#' @details
#' Artificial data was generated following the procedure explained in Box 1 in
#' Sánchez-Pinillos et al. (2023) In particular, the initial state of each community was
#' defined using a hypothetical environmental space with optimal locations for
#' 12 species. Community dynamics were simulated using a general Lotka-Volterra model.
#'
#' State dissimilarities were calculated using the Bray-Curtis metric. Segment and
#' trajectory dissimilarities were calculated using the package 'ecotraj'.
#'
#' @references
#' Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023. Ecological Dynamic
#' Regimes: Identification, characterization, and comparison. *Ecological Monographs*.
#' <doi:10.1002/ecm.1589>
#'
#'
"EDR_data"




