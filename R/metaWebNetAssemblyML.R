#' Default signs of the interaction layers
#'
#' A table with the layers expected by [build_multilayer_assembly()] when the links carry no sign attribute.
#' `sign_from` and `sign_to` are the signs of the effect received by the species at the origin and at the end of
#' each link. The defaults assume that links go from the species with the lower sign to the species with the
#' higher sign: trophic links from the prey (-) to the predator (+), commensalism from the unaffected species (0)
#' to the beneficiary (+), amensalism from the harmed species (-) to the unaffected one (0). Mutualism (+,+) and
#' competition (-,-) are symmetric, and non-trophic antagonism is coded (-,+) like a trophic link.
#'
#' Datasets use different conventions, so check the direction of the links and edit the table (or use a sign
#' attribute per link, see `effect_attr` in [build_multilayer_assembly()]) before simulating.
#'
#' `rule` is `"trophic"` for the layer that defines the meta-web and the prey rule, and `"modifier"` for the
#' layers whose links change the extinction rates.
#'
#' @return a data.frame with columns `layer`, `sign_from`, `sign_to` and `rule`.
#' @seealso [build_multilayer_assembly()], [metaWebNetAssemblyML()]
#' @examples
#' default_layer_signs()
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
#' Translates a multilayer network and the parameters `beta` into the three matrices used by
#' [metaWebNetAssemblyMLCT()]: the trophic meta-web `metaW`, the matrix `M` of effects on the log extinction rate,
#' and the obligate dependence matrix `Dep`. It is called by [metaWebNetAssemblyML()]; use it directly to check how
#' the layers were translated.
#'
#' **Signs of the links.** Each modifier link gives the sign (+1, 0 or -1) of the effect that each of the two
#' species receives from the other. The signs are taken from one of two sources:
#'
#' * a sign attribute per link (`effect_attr`, values `"+"`, `"-"` or numbers): the link `a -> b` is the effect of
#'   `a` on `b` only; an effect of `b` on `a` must be another link `b -> a`. A layer must have the attribute in all
#'   its links or in none.
#' * otherwise, the signs of the layer in `layers`: in a link `a -> b`, `a` receives `sign_from` and `b` receives
#'   `sign_to`.
#'
#' **From signs to `M`.** `M[i, j]` is the effect of the presence of `j` on the log extinction rate of `i`. A positive
#' effect of `j` on `i` adds `-beta["pos"]` to `M[i, j]` and a negative one adds `+beta["neg"]`, both multiplied by
#' `layer_beta[layer]` if given. Pairs linked in several layers add their effects. In symmetric layers
#' (`sign_from == sign_to`) each pair counts once, however many directions are recorded. With `normalize = TRUE` the
#' negative and positive entries of each row are then divided by the number of partners of that sign.
#'
#' **Obligate dependence.** With `obligate_on`, `Dep[i, j]` is `TRUE` if `j` is in `obligate_on`, `i` is not, and
#' `j` has a positive effect on `i` in one of `obligate_layers`.
#'
#' Self-loops are kept in the trophic layer (cannibals are their own prey) and dropped in the modifier layers.
#'
#' @param mg the multilayer network: an `mgraph` (the list of igraph layers returned by `multiweb::readMultiplex()`)
#'   or any list of igraph layers, where the layer id is the list name or, if the list has no names, the
#'   `type_attr` edge attribute; a single igraph with the layer id in the edge attribute `type_attr`; or a
#'   data.frame with columns `from`, `to`, `layer` and, optionally, the sign column `effect_attr`. Vertices need a
#'   `name` attribute.
#' @param layers data.frame with one row per layer: `layer` (the id), `sign_from`, `sign_to` and `rule`
#'   (`"trophic"` or `"modifier"`); see [default_layer_signs()]. Every layer of `mg` must be in it, and exactly one
#'   must have `rule == "trophic"`.
#' @param beta named numeric vector `c(pos = , neg = )`: strength of the positive and of the negative effects on the
#'   log extinction rate. With `normalize = TRUE`, `exp(beta)` is the largest factor by which the partners can
#'   change the extinction rate.
#' @param layer_beta optional named numeric vector, a multiplier of `beta` for each layer (default 1).
#' @param obligate_on character vector of species on which others can depend obligately (e.g. foundation
#'   species), or `NULL` (default) for no obligate dependence.
#' @param obligate_layers layers whose positive effects create obligate dependence; `NULL` (default) uses all
#'   modifier layers.
#' @param normalize if `TRUE` (default) the positive and negative effects on each species are divided by its
#'   number of positive and negative partners, so that the effect depends on the fraction of partners present.
#'   With `FALSE` the effects add up per partner present, and a species with many mutualists can become practically
#'   immortal.
#' @param effect_attr name of the edge attribute (or data.frame column) with the sign of the effect of `from` on
#'   `to`. If it is `NULL` or missing, the layer signs are used.
#' @param type_attr name of the edge attribute with the layer id.
#' @return a list with `metaW` (logical, `metaW[j, i]` is `TRUE` if `j` is a prey of `i`), `M` (numeric), `Dep`
#'   (logical), all with the species as dimnames, and `species`.
#' @seealso [metaWebNetAssemblyML()], [default_layer_signs()]
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

#' Simulation of an assembly process from a multilayer meta-web
#'
#' Continuous time colonization-extinction model of [metaWebNetAssemblyCT()] in which the non-trophic layers
#' change the extinction rates and, optionally, some species depend obligately on others (e.g. on foundation
#' species). Losing a partner with a positive effect (a mutualist, a facilitator) raises the extinction rate of a
#' species; losing a partner with a negative effect (a competitor) lowers it.
#'
#' The extinction rate of species \eqn{i} while species are present is
#' \deqn{e_i(t) = e_i \exp(\beta_- f^-_i(t) - \beta_+ f^+_i(t))}
#' where \eqn{f^+_i(t)} and \eqn{f^-_i(t)} are the fractions of the partners of \eqn{i} with a positive and a
#' negative effect on it that are present at time \eqn{t} (`normalize = TRUE`, the default). Thus
#' \eqn{\exp(\beta_+)} and \eqn{\exp(\beta_-)} are the largest factors by which the partners can lower or raise the
#' extinction rate. Migration and the trophic rule (a predator without prey goes extinct with rate *se*) are as in
#' [metaWebNetAssemblyCT()]. With `obligate_on`, a dependent species goes extinct with rate *se* when none of the
#' species it depends on is present.
#'
#' The parameters `beta` (and `layer_beta`) are translated by [build_multilayer_assembly()] into the matrix `M` of
#' effects on the log extinction rate used by the C++ core [metaWebNetAssemblyMLCT()]. Call
#' [build_multilayer_assembly()] to check how a network was translated.
#'
#' With `beta = c(pos = 0, neg = 0)` and `obligate_on = NULL` the model is [metaWebNetAssemblyCT()] on the trophic layer.
#'
#' @inheritParams build_multilayer_assembly
#' @param m,e,se migration, baseline extinction and secondary extinction rates, one value or one per species,
#'   in the order of the `species` element of the result (the order in which species first appear in `mg`).
#' @param time simulated time, the state is recorded at integer times 1 to time-1.
#' @return the list of [metaWebNetAssemblyMLCT()] (`S`, `L`, `STime`, `A`, `rejected`) with species names in
#'   `STime` and `A`, plus `metaW`, `M`, `Dep` and `species` from [build_multilayer_assembly()].
#' @seealso [build_multilayer_assembly()], [default_layer_signs()], `vignette("multilayer-assembly")`
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
