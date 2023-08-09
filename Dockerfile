# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

FROM genepattern/outsplice:0.4

COPY outsplice_plotjunctions_wrapper.R /build/source/outsplice_plotjunctions_wrapper.R

CMD ["bash"  ]

