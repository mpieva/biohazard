Genotyper
=========

- see which error formula works better: 1/4 or 1/3
- ordering of bases in samtools model: increasing or decrasing quality

- test on simulated and real data (MLE of heterozygosity, separately for
  transitions and transversions, especially on Denise)

Damage-Patterns
===============

The idea is to gather some statistics from BAM files and put the code
for it into biohazard.  Damage-patterns could then be a program that
simply renders the output, and biohazard could contain one that produces
the output as text or JSON.  (JSON is preferred for being structured and
for being understood by flotcharts.)

