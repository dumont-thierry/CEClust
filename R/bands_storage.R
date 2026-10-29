bands_store <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Supply an archive directory.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create archive directory: ", path)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  dir.create(file.path(path, "labels"), showWarnings = FALSE)
  structure(list(path = path, format = 1L), class = "bands_store")
}

bands_lock <- function(path, timeout = 30) {
  lock <- paste0(path, ".lock")
  deadline <- Sys.time() + timeout
  while (!dir.create(lock, showWarnings = FALSE)) {
    if (Sys.time() > deadline) stop("Archive is locked: ", path,
      ". Check that its writer has stopped before removing the lock directory.")
    Sys.sleep(.05)
  }
  lock
}

bands_write_file <- function(value, path) {
  temporary <- tempfile(pattern = ".bands-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = TRUE, version = 3L)
  for (attempt in seq_len(100L)) {
    if (suppressWarnings(file.rename(temporary, path))) return(invisible(path))
    Sys.sleep(.05)
  }
  stop("Cannot atomically publish archive: ", path,
       ". The existing file has not been removed.")
}

bands_save <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- bands_lock(path)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  bands_write_file(value, path)
}

bands_pack_labels <- function(labels) {
  labels <- bands_labels(labels); n <- length(labels); K <- max(labels)
  encoding <- if (K == 1L) "trivial" else if (K == n) "singletons" else
    if (K == 2L) "bits" else if (K <= 255L) "uint8" else
      if (K <= 65535L) "uint16le" else "integer"
  values <- switch(encoding,
    trivial = NULL, singletons = NULL,
    bits = packBits(as.raw(c(labels - 1L, rep(0L, (-n) %% 8L))), type = "raw"),
    uint8 = as.raw(labels),
    uint16le = as.raw(as.vector(rbind(bitwAnd(labels, 255L), bitwShiftR(labels, 8L)))),
    integer = as.integer(labels))
  list(format = 1L, n = n, K = K, encoding = encoding, values = values)
}

bands_unpack_labels <- function(value) {
  if (!identical(value$format, 1L)) stop("Unknown allocation archive format.")
  labels <- switch(value$encoding,
    trivial = rep(1L, value$n), singletons = seq_len(value$n),
    bits = as.integer(rawToBits(value$values))[seq_len(value$n)] + 1L,
    uint8 = as.integer(value$values),
    uint16le = {
      bytes <- matrix(as.integer(value$values), nrow = 2L)
      as.integer(bytes[1L, ] + 256L * bytes[2L, ])
    },
    integer = value$values,
    stop("Unknown allocation encoding."))
  if (length(labels) != value$n || anyNA(labels) || any(labels < 1L) ||
      max(labels) != value$K) stop("Corrupt allocation archive.")
  as.integer(labels)
}

bands_candidate <- function(X, labels, store, origin = "", data_key = NULL) {
  X <- bands_data(X)
  if (length(labels) != nrow(X) || anyNA(labels)) stop("Invalid candidate allocation.")
  labels <- bands_labels(labels)
  if (is.null(data_key)) data_key <- digest::digest(X, algo = "sha256")
  if (length(data_key) != 1L || !is.character(data_key) || is.na(data_key) || !nzchar(data_key))
    stop("data_key must identify the complete ordered training sample.")
  label_key <- digest::digest(labels, algo = "xxhash64")
  key <- digest::digest(list(data_key, label_key), algo = "xxhash64")
  relative <- file.path("labels", paste0(key, ".rds")); path <- file.path(store$path, relative)
  lock <- bands_lock(path)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  packed <- NULL
  if (file.exists(path)) {
    packed <- readRDS(path)
    if (!identical(packed$data_key, data_key) || !identical(packed$label_key, label_key) ||
        !identical(bands_unpack_labels(packed), labels)) stop("Candidate hash collision: ", key)
  }
  if (is.null(packed) || !identical(packed$moments_format, 1L)) {
    packed <- bands_pack_labels(labels)
    packed$data_key <- data_key; packed$label_key <- label_key
    packed$rule <- bands_rule(X, labels)
    packed$stats <- bands_rule_stats(packed$rule)
    packed$moments_format <- 1L
    bands_write_file(packed, path)
  }
  if (!identical(packed$rule$n, nrow(X)) || !identical(packed$rule$d, ncol(X)) ||
      nrow(packed$stats) != max(labels)) stop("Candidate moments archive mismatch: ", key)
  list(key = key, id = key, data_key = data_key, label_key = label_key,
       n = nrow(X), d = ncol(X), K = max(labels), stats = packed$stats,
       rule = packed$rule, label_file = relative, origin = origin)
}

bands_load_labels <- function(candidate, store) {
  path <- file.path(store$path, candidate$label_file)
  # Published allocation files are immutable: atomic replacement only adds
  # equivalent cached moments. Readers do not acquire a writing lock.
  value <- bands_read_file(path)
  if (!identical(value$data_key, candidate$data_key) ||
      !identical(value$label_key, candidate$label_key)) stop("Candidate archive identity mismatch.")
  labels <- bands_unpack_labels(value)
  if (!identical(digest::digest(labels, algo = "xxhash64"), candidate$label_key))
    stop("Candidate allocation checksum mismatch.")
  labels
}

bands_read_file <- function(path) {
  for (attempt in seq_len(20L)) {
    value <- try(readRDS(path), silent = TRUE)
    if (!inherits(value, "try-error")) return(value)
    if (!file.exists(path)) break
    Sys.sleep(.05)
  }
  stop("Cannot read archive: ", path, "\n", as.character(value))
}

bands_cache <- function(max_bytes = 64 * 1024^2) {
  if (length(max_bytes) != 1L || !is.finite(max_bytes) || max_bytes < 0)
    stop("Cache budget must be a finite nonnegative number of bytes.")
  cache <- new.env(parent = emptyenv())
  cache$values <- new.env(hash = TRUE, parent = emptyenv())
  cache$sizes <- new.env(hash = TRUE, parent = emptyenv())
  cache$order <- character(); cache$bytes <- 0; cache$max_bytes <- max_bytes
  cache$hits <- 0L; cache$misses <- 0L; cache$evictions <- 0L
  cache
}

bands_cache_get <- function(cache, key) {
  if (is.null(cache)) return(NULL)
  if (!exists(key, cache$values, inherits = FALSE)) {
    cache$misses <- cache$misses + 1L
    return(NULL)
  }
  cache$hits <- cache$hits + 1L
  cache$order <- c(cache$order[cache$order != key], key)
  get(key, cache$values, inherits = FALSE)
}

bands_cache_put <- function(cache, key, value) {
  if (is.null(cache)) return(invisible(FALSE))
  size <- as.numeric(utils::object.size(value)) + as.numeric(utils::object.size(key)) + 256
  if (size > cache$max_bytes) return(invisible(FALSE))
  if (exists(key, cache$values, inherits = FALSE)) bands_cache_remove(cache, key)
  while (cache$bytes + size > cache$max_bytes && length(cache$order)) {
    bands_cache_remove(cache, cache$order[1L]); cache$evictions <- cache$evictions + 1L
  }
  assign(key, value, cache$values); assign(key, size, cache$sizes)
  cache$order <- c(cache$order, key); cache$bytes <- cache$bytes + size
  invisible(TRUE)
}

bands_cache_remove <- function(cache, key) {
  cache$bytes <- cache$bytes - get(key, cache$sizes, inherits = FALSE)
  rm(list = key, envir = cache$values); rm(list = key, envir = cache$sizes)
  cache$order <- cache$order[cache$order != key]
  invisible(NULL)
}
