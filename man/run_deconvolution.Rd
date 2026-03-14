% Generated manually for P9
\name{run_deconvolution}
\alias{run_deconvolution}
\alias{run_aegis_full}
\alias{run_spotlight}
\alias{run_card}
\alias{run_rctd}
\alias{run_cell2location}
\alias{run_destvi}
\alias{run_tangram}
\alias{run_stereoscope}
\title{Unified deconvolution orchestration and method runners}
\usage{
run_deconvolution(
  seu,
  reference,
  methods = c("SPOTlight", "RCTD", "CARD"),
  sample_id = NULL,
  normalize = TRUE,
  strict = TRUE,
  use_python = TRUE,
  runner_overrides = NULL,
  runner_args = NULL,
  ...
)

run_aegis_full(
  seu,
  reference,
  methods = c("SPOTlight", "RCTD", "CARD"),
  markers = NULL,
  strict = TRUE,
  use_python = TRUE,
  normalize = TRUE,
  ...
)

run_spotlight(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_card(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_rctd(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_cell2location(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_destvi(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_tangram(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_stereoscope(seu, reference, normalize = TRUE, strict = TRUE, ...)
}
\arguments{
\item{seu}{A Seurat spatial object.}

\item{reference}{Reference input used by runnable methods.}

\item{methods}{Character vector of method names.}

\item{sample_id}{Optional sample identifier recorded in deconvolution output metadata.}

\item{normalize}{Logical; row-normalize outputs to proportions when possible.}

\item{strict}{Logical; strict error behavior for unsupported/missing dependencies or failed execution.}

\item{use_python}{Logical; allow Python-backed methods via \code{reticulate}.}

\item{runner_overrides}{Optional named list of runner functions keyed by method, mainly for advanced usage/testing.}

\item{runner_args}{Optional named list of per-method argument lists passed to the selected runner.}

\item{markers}{Optional marker list used by \code{run_aegis_full()}.}

\item{...}{Additional arguments forwarded to method runners.}
}
\value{
\code{run_deconvolution()} returns a list with:\cr
\code{seu}, \code{deconv}, \code{methods_run}, \code{methods_skipped}, and \code{messages}.\cr
\code{run_aegis_full()} returns an \code{aegis} object.\cr
Method-specific runners return spot-by-celltype numeric matrices.
}
\description{
\code{run_deconvolution()} dispatches requested methods using the AEGIS method registry,
standardizes outputs to the shared AEGIS deconvolution contract, and returns
results ready for \code{as_aegis()} / \code{run_aegis()}.

\code{run_aegis_full()} combines one-shot deconvolution dispatch with the existing
AEGIS downstream pipeline.

Method-specific runner wrappers are included for directly runnable R-native methods
and optional Python-backed methods.
}
