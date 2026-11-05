# Load named workflow modules, set collation and provide seed and log helpers.
# Rules declare their consumed modules as file inputs. Sourcing uses base R;
# parse_script_arguments() needs optparse only when called.

PIPELINE_SEED <- 1337L

#' Set LC_COLLATE=C so sorting and clonotype suffixes agree across runtimes.
#'
#' Stops if the locale cannot be set. Returns nothing.
set_pipeline_collation <- function() {
  if (!nzchar(Sys.setlocale("LC_COLLATE", "C"))) {
    stop("Cannot set LC_COLLATE=C; sort order would depend on the runtime image")
  }
  invisible(NULL)
}

#' Find the shared R directory from a directory within the workflow checkout.
#'
#' from defaults to the running script's directory, read from Rscript's --file
#' argument. Callers without that argument must supply from. Returns an absolute
#' path, or stops when the workflow root cannot be found.
pipeline_r_directory <- function(from = NULL) {
  if (is.null(from)) {
    scripts <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(scripts) != 1L) {
      stop("Supply from when no running script path is available")
    }
    from <- dirname(normalizePath(sub("^--file=", "", scripts[[1L]])))
  }
  candidate <- normalizePath(from, mustWork = TRUE)
  repeat {
    modules <- file.path(candidate, "workflow", "lib")
    if (file.exists(file.path(candidate, "Snakefile")) &&
        file.exists(file.path(modules, "runtime.R"))) {
      return(normalizePath(modules))
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop("Cannot find the workflow root from ", from)
}

#' Load the named shared modules.
#'
#' `modules` are base names without the .R suffix, in the order they should be
#' sourced; a module that calls into another must come after it. `from` is a
#' directory the workflow root is reachable from, or NULL for a stage script,
#' which recovers its own location. Returns nothing. Side effect: defines each
#' module's functions and constants in the global environment.
load_pipeline_modules <- function(modules, from = NULL) {
  set_pipeline_collation()
  r_directory <- pipeline_r_directory(from)
  for (module in modules) {
    path <- file.path(r_directory, paste0(module, ".R"))
    if (!file.exists(path)) {
      stop("No shared R module named '", module, "'. Looked for ", path)
    }
    source(path)
  }
  invisible(NULL)
}

# The shell redirects stdout and stderr to the rule's log file. These helpers
# write messages to stderr without opening a file.
LOG_TIMESTAMP_FORMAT <- "%Y-%m-%d %H:%M:%S"

#' Write one timestamped line at the given level.
#'
#' `format_string` and `...` are passed to sprintf, so a message with no
#' substitutions can be given on its own. Writes to stderr through message().
log_line <- function(level, format_string, ...) {
  message(sprintf("%s [%s] %s", format(Sys.time(), LOG_TIMESTAMP_FORMAT), level,
                  sprintf(format_string, ...)))
}

#' Write one line about the stage's progress.
log_info <- function(format_string, ...) log_line("INFO", format_string, ...)

#' Write a recoverable condition to stderr without signalling an R warning.
#'
#' Use warning() when condition handlers need the warning, or stop() when the
#' output cannot be trusted. start_stage_log() makes warnings print immediately.
log_warn <- function(format_string, ...) log_line("WARN", format_string, ...)

#' Log the stage's start and record its start time.
#'
#' `stage_name` is the entry-point script's name without its extension.
#'
#' Returns a small list to hand back to finish_stage_log(). Side effects: writes
#' the opening line and makes warnings print where they happen so the log stays
#' in order.
start_stage_log <- function(stage_name) {
  options(warn = 1)
  log_info("Started %s", stage_name)
  list(stage_name = stage_name, started_at = Sys.time())
}

#' Log the stage's elapsed time.
#'
#' `stage_log` is the value from start_stage_log(). Returns nothing.
finish_stage_log <- function(stage_log) {
  elapsed_seconds <- as.numeric(
    difftime(Sys.time(), stage_log$started_at, units = "secs"))
  log_info("Completed %s in %.1f s", stage_log$stage_name, elapsed_seconds)
  invisible(NULL)
}

#' Parse this script's command-line arguments.
#'
#' `option_list` holds optparse::make_option() definitions. An option given no
#' default is required, and a run without it stops after printing the usage.
#'
#' Returns the named list of resolved values and logs them. The log does not
#' quote shell arguments or include the script name.
parse_script_arguments <- function(option_list) {
  parser <- optparse::OptionParser(option_list = option_list)
  script_options <- optparse::parse_args(parser)

  required <- vapply(option_list, function(option) option@dest, character(1))[
    vapply(option_list, function(option) is.null(option@default), logical(1))
  ]
  absent <- setdiff(required, names(script_options))
  if (length(absent)) {
    optparse::print_help(parser)
    stop("Missing required argument(s): --", paste(absent, collapse = ", --"))
  }

  supplied <- script_options[setdiff(names(script_options), "help")]
  log_info("Arguments: %s", paste(sprintf("--%s=%s", names(supplied),
                                          vapply(supplied, format, character(1))),
                                  collapse = " "))
  script_options
}

#' Split a comma-separated argument into paths or identifiers in input order.
#'
#' Removes surrounding whitespace and empty components. Stops if none remain.
#' Commas cannot occur within a path or identifier in this encoding.
split_path_list <- function(value) {
  paths <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  paths <- paths[nzchar(paths)]
  if (!length(paths)) stop("Expected at least one path in a comma-separated list")
  paths
}

#' Parse comma-separated NAME=value pairs into a named numeric vector.
#'
#' Keeps input order and ignores blank components. An empty argument returns an
#' empty named vector. Stops on malformed pairs, blank or repeated names and
#' non-numeric or missing values.
parse_named_values <- function(value) {
  components <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  components <- components[nzchar(components)]
  if (!length(components)) return(setNames(numeric(0), character(0)))
  pairs <- regmatches(components, regexec("^([^=]+)=([^=]*)$", components))
  malformed <- vapply(pairs, function(pair) length(pair) != 3L, logical(1))
  if (any(malformed)) {
    stop("Expected NAME=value components, got: ",
         paste(components[malformed], collapse = ", "))
  }
  names <- trimws(vapply(pairs, `[[`, character(1), 2L))
  values <- suppressWarnings(as.numeric(vapply(pairs, `[[`, character(1), 3L)))
  if (any(!nzchar(names)) || anyNA(values)) {
    stop("Expected NAME=number components, got: ", paste(components, collapse = ", "))
  }
  if (anyDuplicated(names)) {
    stop("Repeated name in named values: ", paste(names[duplicated(names)], collapse = ", "))
  }
  setNames(values, names)
}
