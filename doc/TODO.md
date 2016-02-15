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

Options for damage-patterns:

  -a       --adapter                        Require B adaptor even for 5' graphs
        Amounts to "short reads only" and is generally uninteresting.  

  -A       --no-adapter                     Require absence of adapter
        Amounts to "long reads only" and is generally uninteresting.
        
  -5       --five                           Treat 5' end only
        Useless.

  -b NUM   --min-length=NUM, --bin=NUM      Only consider length bin starting at NUM
        Not generally useful.

  -B NUM   --max-length=NUM, --bin-end=NUM  Only consider length bin ending at NUM
        Not generally useful.

  -I FILE  --filter=FILE                    Use only ids from file FILE
        Not generally useful.

  -X NUM   --chop=NUM                       Chop off NUM nt at far end
        Silly, nobody uses it.

  -x NUM   --xrange=NUM                     Set range of X axis to [0..NUM]
        Silly, nobody uses it.

  -C NUM   --context=NUM                    Extract NUM nt context from reference
        Doesn't need to be configurable, a default of 10 is bloody enough.

  -G FILE  --genome=FILE                    Genome in 2bit format to get context/ref from
        Needed for context, needed in absence of MD field.

           --no-cpg                         Disregard CpG motifs
        Worth thinking about.

           --bylength                       Additionally create deamination-by-length-bin graphs
        Might be useful and reasonably easy to support.
    
           --leeHom                         Pretend everything is adapter-trimmed
        Needed at this time, but not generally useful.

