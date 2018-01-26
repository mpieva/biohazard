-- | Umbrella module for most of what's under 'Bio.Bam'.

module Bio.Bam (
    module Bio.Bam.Fastq,
    module Bio.Bam.Filter,
    module Bio.Bam.Header,
    module Bio.Bam.Index,
    module Bio.Bam.Reader,
    module Bio.Bam.Rec,
    module Bio.Bam.Trim,
    module Bio.Bam.Writer,
    module Bio.Iteratee
               ) where

import Bio.Bam.Fastq
import Bio.Bam.Filter
import Bio.Bam.Header
import Bio.Bam.Index
import Bio.Bam.Reader
import Bio.Bam.Rec
import Bio.Bam.Trim
import Bio.Bam.Writer
import Bio.Iteratee

