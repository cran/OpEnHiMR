\name{pred_hmc}
\alias{pred_hmc}
\title{Prediction of Histone Modification for Multiclass Data}
\usage{
pred_hmc(fasta_file_path)
}
\arguments{
\item{fasta_file_path}{Sequence file path (.fasta format)}
}
\value{
Modications: sequences with their modifications (H3K27me3, H3K9ac and H3K4me3) or no modification.
}
\description{
Prediction of H3K27me3, H3K9ac and H3K4me3 modification in rice.
}
\examples{
\donttest{
example_path <- system.file("exdata/test.fasta", package = "OpEnHiMR")
pred <- pred_hmc(fasta_file_path = example_path)
}
}
\references{
Yin, Q., Wu, M., Liu, Q. et al. DeepHistone: a deep learning approach to predicting histone modifications. BMC Genomics 20 (Suppl 2), 193 (2019).
}
