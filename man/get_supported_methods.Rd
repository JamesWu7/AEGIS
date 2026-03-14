% Generated manually for P9
\name{get_supported_methods}
\alias{get_supported_methods}
\title{List supported deconvolution methods and execution modes}
\usage{
get_supported_methods()
}
\value{
A \code{data.frame} with method capabilities, support mode, dependency type,
and adapter/runner mapping.
}
\description{
Returns the AEGIS method registry used by \code{run_deconvolution()}.
The registry explicitly separates direct R execution, optional Python execution,
and import-only methods.
}
