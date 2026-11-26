#' Default signs of the interaction layers
#'
#' Each layer is coded by the sign of the effect on the species at the origin (`from`) and at the end (`to`)
#' of its links. The convention is that links go from the species with the lower sign to the species with the
#' higher sign: trophic links from the prey (-) to the predator (+), commensalism from the unaffected species (0)
#' to the beneficiary (+), amensalism from the harmed species (-) to the unaffected one (0). Mutualism and
#' competition are symmetric.
#'
#' `rule` is `"trophic"` for the layer used for the prey rule (a predator without prey goes extinct with rate *se*)
#' and `"modifier"` for the layers whose links change the extinction rates.
#'
#' Non-trophic antagonism is coded as (-,+), as a trophic link; change it if that is not what the data mean.
#'
#' @return a data.frame with columns `layer`, `sign_from`, `sign_to` and `rule`.
#' @seealso [metaWebNetAssemblyML()]
#' @export
default_layer_signs <- function() {
  data.frame(layer     = c("Trophic", "Antagonism", "Mutualism", "Competition", "Commensalism", "Amensalism"),
             sign_from = c(-1, -1,  1, -1, 0, -1),
             sign_to   = c( 1,  1,  1, -1, 1,  0),
             rule      = c("trophic", "modifier", "modifier", "modifier", "modifier", "modifier"),
             stringsAsFactors = FALSE)
}

# Edge list (from, to, layer, effect) from an mgraph / named list of igraph layers, a single igraph with a type
# attribute, or a data.frame with columns from, to and layer.
.ml_edges <- function(mg, type_attr, effect_attr = NULL) {
  ig_edges <- function(g, layer = NULL) {
    if (!requireNamespace("igraph", quietly = TRUE)) stop("package igraph is needed for igraph input")
    if (is.null(igraph::V(g)$name)) stop("the vertices of every layer need a name attribute")
    d <- igraph::as_data_frame(g, what = "edges")
    if (is.null(layer)) {
      if (!type_attr %in% names(d)) stop("edge attribute '", type_attr, "' not found")
      layer <- d[[type_attr]]
    }
    eff <- if (!is.null(effect_attr) && effect_attr %in% names(d)) as.character(d[[effect_attr]]) else NA_character_
    list(e = data.frame(from = d$from, to = d$to, layer = as.character(layer), effect = eff, stringsAsFactors = FALSE),
         sp = igraph::V(g)$name)
  }
  if (is.data.frame(mg)) {
    stopifnot(all(c("from", "to", "layer") %in% names(mg)))
    e <- data.frame(from = as.character(mg$from), to = as.character(mg$to), layer = as.character(mg$layer),
                    effect = if (!is.null(effect_attr) && effect_attr %in% names(mg)) as.character(mg[[effect_attr]]) else NA_character_,
                    stringsAsFactors = FALSE)
    return(list(e = e, sp = unique(c(e$from, e$to))))
  }
  if (inherits(mg, "igraph")) return(ig_edges(mg))
  if (is.list(mg)) {
    nm <- names(mg)
    parts <- lapply(seq_along(mg), function(k) {
      g <- mg[[k]]
      lay <- if (!is.null(nm) && nzchar(nm[k])) nm[k] else NULL       # list name is the layer id
      if (is.null(lay)) {                                             # otherwise the type attribute
        ty <- unique(igraph::edge_attr(g, type_attr))
        if (length(ty) != 1) stop("layer ", k, " has no name and no single '", type_attr, "' value")
        lay <- ty
      }
      ig_edges(g, rep(lay, igraph::ecount(g)))
    })
    return(list(e = do.call(rbind, lapply(parts, `[[`, "e")), sp = unique(unlist(lapply(parts, `[[`, "sp")))))
  }
  stop("mg must be an mgraph or list of igraph layers, an igraph, or a data.frame (from, to, layer)")
}

#' Matrices of the multilayer assembly model
#'
#' Builds the trophic matrix, the extinction modifier matrix and the obligate dependence matrix used by
#' [metaWebNetAssemblyMLCT()] from a multilayer network. Called by [metaWebNetAssemblyML()]; useful by itself to
#' check how the layers were translated.
#'
#' If the links have a sign attribute (`effect_attr`, values `"+"`/`"-"` or numbers), each link `a -> b` is read
#' as the effect of `a` on `b` and only adds `-beta_pos` (positive) or `+beta_neg` (negative) to `M[b, a]`; the effect
#' of `b` on `a`, if any, must be another link `b -> a`. This is how the mussel-fucus networks are coded, and then
#' the signs in `layers` are not used for the modifier layers. Otherwise the layer signs are used:
#' for a modifier layer with signs (`sign_from`, `sign_to`), a link `a -> b` adds `-beta_l * sign_to` to `M[b, a]`
#' and `-beta_l * sign_from` to `M[a, b]`, where `beta_l` is `beta["pos"]` for a positive sign and `beta["neg"]`
#' for a negative one, times `layer_beta[layer]`. Thus the presence of a species with a positive effect lowers the
#' extinction rate and one with a negative effect raises it. Pairs linked in several layers add their effects.
#' In symmetric layers each pair is counted once whatever the number of directions recorded.
#'
#' @param mg the multilayer network: an `mgraph` or a list of igraph layers (layer id = list name, or the
#'   `type_attr` edge attribute if the list has no names), a single igraph with the layer in the edge attribute
#'   `type_attr`, or a data.frame with columns `from`, `to` and `layer`.
#' @param layers data.frame with the layer signs and rules, see [default_layer_signs()]. Every layer present in
#'   `mg` must be in it, and exactly one layer must have `rule == "trophic"`.
#' @param beta numeric vector `c(pos = , neg = )`, strength of positive and negative effects on the log extinction rate.
#' @param layer_beta optional named numeric vector, multiplies `beta` for each layer (default 1).
#' @param obligate_on character vector of species (e.g. the foundation species). A species that receives a
#'   positive effect from any of them in the layers `obligate_layers` depends obligately on them (and is not one
#'   of them). `NULL` for no obligate dependence.
#' @param obligate_layers layers whose positive effects create obligate dependence; `NULL` (default) uses all
#'   modifier layers.
#' @param normalize if `TRUE` (default) the positive and the negative effects on each species are divided by
#'   its number of positive and negative partners in the meta-web, so the extinction rate is
#'   \eqn{e_i \exp(\beta_- f^-_i - \beta_+ f^+_i)}, where \eqn{f} are the fractions of partners present, and
#'   \eqn{\exp(\beta)} is the largest possible change. With `FALSE` the effects add up per partner present, and a
#'   species with many mutualists can become practically immortal.
#' @param effect_attr name of the edge attribute (or data.frame column) with the sign of the effect of `from` on
#'   `to`; `NULL` or a missing attribute to use the layer signs.
#' @param type_attr name of the edge attribute with the layer id.
#' @return a list with `metaW`, `M`, `Dep` (species as dimnames) and `species`.
#' @export
build_multilayer_assembly <- function(mg, layers = default_layer_signs(), beta = c(pos = 0, neg = 0),
                                      layer_beta = NULL, obligate_on = NULL, type_attr = "type",
                                      effect_attr = "effect", obligate_layers = NULL, normalize = TRUE) {
  x <- .ml_edges(mg, type_attr, effect_attr)
  E <- x$e; sp <- x$sp; S <- length(sp)
  miss <- setdiff(unique(E$layer), layers$layer)
  if (length(miss)) stop("layers not described in 'layers': ", paste(miss, collapse = ", "))
  troph <- layers$layer[layers$rule == "trophic"]
  if (length(troph) != 1) stop("exactly one layer must have rule 'trophic'")
  if (!all(c("pos", "neg") %in% names(beta))) stop("beta must be c(pos = , neg = )")

  metaW <- matrix(FALSE, S, S, dimnames = list(sp, sp))
  Et <- E[E$layer == troph, ]
  metaW[cbind(Et$from, Et$to)] <- TRUE                   # cannibalism loops are kept

  M <- matrix(0, S, S, dimnames = list(sp, sp))
  P <- matrix(FALSE, S, S, dimnames = list(sp, sp))      # P[i, j]: j has a positive effect on i
  coef <- function(s, lay) {
    lb <- if (!is.null(layer_beta) && lay %in% names(layer_beta)) layer_beta[[lay]] else 1
    ifelse(s > 0, -beta[["pos"]], ifelse(s < 0, beta[["neg"]], 0)) * lb
  }
  Pl <- list()                                           # positive effects by layer, for obligate dependence
  for (lay in setdiff(unique(E$layer), troph)) {
    L <- layers[layers$layer == lay, ]
    if (L$rule != "modifier") stop("rule of layer ", lay, " must be 'trophic' or 'modifier'")
    El <- E[E$layer == lay & E$from != E$to, ]
    Pk <- matrix(FALSE, S, S, dimnames = list(sp, sp))
    if (nrow(El) && all(!is.na(El$effect))) {            # signed links: effect of from on to
      El <- El[!duplicated(paste(El$from, El$to, sep = "\r")), ]
      sg <- ifelse(El$effect %in% c("+", "1", "+1"), 1, ifelse(El$effect %in% c("-", "-1"), -1, suppressWarnings(as.numeric(El$effect))))
      if (anyNA(sg)) stop("unknown values in the effect attribute of layer ", lay)
      for (k in seq_len(nrow(El))) {
        M[El$to[k], El$from[k]] <- M[El$to[k], El$from[k]] + coef(sign(sg[k]), lay)
        if (sg[k] > 0) Pk[El$to[k], El$from[k]] <- TRUE
      }
      Pl[[lay]] <- Pk
      next
    }
    if (nrow(El) && any(!is.na(El$effect))) stop("layer ", lay, " has the effect attribute only in some links")
    if (L$sign_from == L$sign_to) {                      # symmetric: one link per pair
      key <- paste(pmin(El$from, El$to), pmax(El$from, El$to), sep = "\r")
      El <- El[!duplicated(key), ]
    } else El <- El[!duplicated(paste(El$from, El$to, sep = "\r")), ]
    if (!nrow(El)) next
    ab <- cbind(El$to, El$from); ba <- cbind(El$from, El$to)
    M[ab] <- M[ab] + coef(L$sign_to, lay)                 # effect of 'from' on 'to'
    M[ba] <- M[ba] + coef(L$sign_from, lay)               # effect of 'to' on 'from'
    if (L$sign_to > 0)   Pk[ab] <- TRUE
    if (L$sign_from > 0) Pk[ba] <- TRUE
    Pl[[lay]] <- Pk
  }
  if (normalize) {                                       # effect = beta * fraction of partners present
    Mp <- pmin(M, 0); Mn <- pmax(M, 0)
    np <- rowSums(Mp < 0); nn <- rowSums(Mn > 0)
    Mp <- Mp / pmax(np, 1); Mn <- Mn / pmax(nn, 1)
    M <- Mp + Mn
  }
  ol <- if (is.null(obligate_layers)) names(Pl) else intersect(obligate_layers, names(Pl))
  for (lay in ol) P <- P | Pl[[lay]]

  Dep <- matrix(FALSE, S, S, dimnames = list(sp, sp))
  if (!is.null(obligate_on)) {
    bad <- setdiff(obligate_on, sp)
    if (length(bad)) stop("species in obligate_on not found: ", paste(bad, collapse = ", "))
    Dep[, obligate_on] <- P[, obligate_on, drop = FALSE]
    Dep[obligate_on, ] <- FALSE
  }
  list(metaW = metaW, M = M, Dep = Dep, species = sp)
}

#' Simulation of an Assembly process from a multilayer Meta-web
#'
#' Continuous time colonization-extinction model of [metaWebNetAssemblyCT()] in which the non-trophic layers
#' change the extinction rates and, optionally, some species depend obligately on others (e.g. on the foundation
#' species). Losing a mutualist or a facilitator raises the extinction rate of a species, losing a competitor
#' lowers it. See [build_multilayer_assembly()] for how the layers are translated and [metaWebNetAssemblyMLCT()]
#' for the dynamics.
#'
#' With `beta = c(pos = 0, neg = 0)` and `obligate_on = NULL` the model is [metaWebNetAssemblyCT()] on the trophic layer.
#'
#' @inheritParams build_multilayer_assembly
#' @param m,e,se migration, baseline extinction and secondary extinction rates, one value or one per species
#'   (in the order of the `species` element of the result).
#' @param time simulated time, the state is recorded at integer times 1 to time-1.
#' @return the list of [metaWebNetAssemblyMLCT()] with species names in `STime` and `A`, plus the matrices
#'   `metaW`, `M` and `Dep` and the vector `species`.
#' @examples
#' ed <- data.frame(from  = c("alga", "grazer", "alga",       "grazer2", "grazer"),
#'                  to    = c("grazer", "fish",  "grazer2",   "fish",    "grazer2"),
#'                  layer = c("Trophic", "Trophic", "Trophic", "Trophic", "Competition"))
#' r <- metaWebNetAssemblyML(ed, beta = c(pos = 1, neg = 1), m = 0.1, e = 0.05, se = 0.5, time = 100)
#' r$STime[, 99]
#' @export
metaWebNetAssemblyML <- function(mg, layers = default_layer_signs(), beta = c(pos = 0, neg = 0),
                                 m, e, se, time, layer_beta = NULL, obligate_on = NULL, type_attr = "type",
                                 effect_attr = "effect", obligate_layers = NULL, normalize = TRUE) {
  b <- build_multilayer_assembly(mg, layers, beta, layer_beta, obligate_on, type_attr, effect_attr,
                                 obligate_layers, normalize)
  S <- length(b$species)
  rep_s <- function(v, nm) { if (length(v) == 1) rep(v, S) else if (length(v) == S) v else
    stop(nm, " must have length 1 or one value per species") }
  res <- metaWebNetAssemblyMLCT(b$metaW, b$M, b$Dep, rep_s(m, "m"), rep_s(e, "e"), rep_s(se, "se"), time)
  rownames(res$STime) <- b$species
  dimnames(res$A) <- list(b$species, b$species)
  c(res, b)
}
