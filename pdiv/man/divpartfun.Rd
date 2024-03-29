% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addpart.R
\name{addpart}
\alias{addpart}
\alias{tripart}
\title{Additive and Tripartite Partitioning}
\usage{
addpart(
  depmix,
  data,
  name.CE = NULL,
  name.SE = NULL,
  excl.zeroes = FALSE,
  fractions = NULL,
  groups = NULL
)

tripart(
  depmix,
  data,
  name.TIC = NULL,
  name.TDC = NULL,
  name.DOM = NULL,
  excl.zeroes = FALSE,
  groups = NULL
)
}
\arguments{
\item{depmix}{A 'model formula' describing the structure of the
data, in the \code{form y ~ mixture/species + plot}}

\item{data}{Data frame describing the composition of plant
communities, on a per-species basis}

\item{name.CE}{Name of column holding the calculated complemetarity
effect. Defaults to \code{CE.y}, where \code{y} is replaced by
the name of the dependent variable (usually the biomass
measure).}

\item{name.SE}{Name of column holding the calculated selection
effect. Defaults to \code{SE.y}, where \code{y} is replaced by
the name of the dependent variable specified.}

\item{excl.zeroes}{logical specifying whether species with unknown
or zero biomass are excluded from calculations (default \code{FALSE}).
Data that are NA are also excluded.}

\item{fractions}{fraction at which species were established in the
community. Normally, this is the inverse of species richness.}

\item{groups}{The calculations can be performed by groups specified
with a right hand side only model formula.}

\item{name.TIC}{Name of column holding the calculated
trait-independent complemetarity effect. Defaults to
\code{TIC.y}, where \code{y} is replaced by the name of the
dependent variable (usually the biomass measure). Numerically
this term is identical to the complementarity effect of the
additive partitioning scheme.}

\item{name.TDC}{Name of column holding the calculated
trait-dependent complementarity effect. Defaults to
\code{TDC.y}, where \code{y} is replaced by the name of the
dependent variable specified.}

\item{name.DOM}{Name of column holding the calculated selection
dominance effect. Defaults to \code{DOM.y}, where \code{y} is
replaced by the name of the dependent variable specified.}
}
\value{
A data frame containing columns for plant species mixture,
    plot, and the computed complementarity and selection effects.
}
\description{
Given a data set that specifies the composition of plant
communities on a per species basis, these functions compute Loreau and
Hector's additive partitioning, or Fox' tripartite partitioning.
}
\details{
The input data frame needs to contain columns defining the units
for which complementarity and selection effects are to be
calculated (typically field plots), the community composition of
these plots, the species contained in these plots, plus their
contribution to the community-level metric that is analysed
(typically biomass). These columns are specified using a model
formula of the form \code{form y ~ mixture/species + plot}. Note
that plot and mixture can be identical, in which case one has to
write \code{form y ~ mixture/species + mixture}.

Species composition is a character string constructed so that the
presence of individual species can be tested by searching for the
occurence of the species code in the composition string. For
example, species \code{A} will be found in composition \code{ABCD}
but not in \code{EFGH}. The composition codes can contain
additional characters as long as the above rule is met,
e.g. \code{A|B|C|D} would work as well.

Note: Monocultures are identified by testing if their composition
codes do not contain any of the other composition codes. It
therefore is crucial to include the monocultures in the data set,
even if their biomass is zero or not known. Otherwise, mixtures may
be misidentified as monocultures (because there are no
corresponding lines for their component monocultures). Example:
composition 'AB' will be treated as monoculture if no composition
'A' or 'B' are present in the data set.

Monoculture biomass (the reference for the underlying relative
yield calculation) is calculated as mean of all the monoculture
plots, excluding \code{NA}s.

Relative yields, on which the additive partitioning scheme bases,
cannot be calculated for communities that contain species with zero
or unknown monoculture biomass. For such plots, the complementarity
and selection effects will be set to \code{NA}, unless
\code{excl.zeroes=TRUE} is specified. Then, these species will be
removed from the species set before complementarity and selection
effects are calculated. The underlying relative yield data are
still calculated using the original species numbers; only then
are the calculated CE plus SE equal to the net effect calculated
as mixture yield minus average monoculture yield (incl. the failing
monoculture with zero yield).

For applications in which species were not initially established in
the same proportions, the respective \code{fractions} can be
provided. In this case, an adjustment of the biomass values is made
prior to the calculations, i.e. the additive partitioning is
performed as if all species were established at the same initial
density. Note that then the net diversity effect CE + SE is no
longer equal to the overyielding of the particular plot with the
non-equal initial proportions!


These partitioning schemes can be performed by groups
(e.g. different experimental treatments or dates). These groups are
specified as right-hand side only formula. If the code runs on a
platform other than Windows and the library \code{parallel} is
installed, the processing of the different groups is parallelized
on different processor cores.
}
\examples{

data(divpart)
# additive partitioning
addpart(m ~ comp / sp + plot, data = divpart, groups = ~group)
# tripartite partitioning
tripart(m ~ comp / sp + plot, data = divpart, groups = ~group)

}
\references{
Loreau M \& Hector A (2001) Partitioning selection and
    complementarity in biodiversity experiments. Nature 412, 72-76

Fox (2005) Interpreting the selection effect of
    biodiversity on ecosystem functioning. Ecology Letters 8,
    846--856
}
\author{
Pascal Niklaus \email{pascal.niklaus@ieu.uzh.ch}
}
