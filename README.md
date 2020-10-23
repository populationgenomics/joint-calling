# Variant QC

Variant and sample QC, using [gnomad_qc tools](https://github.com/broadinstitute/gnomad_qc)

The gnomad_qc README states:

 > The scripts make reference to gnomAD-related metadata files (not public) and may perform procedures that are not strictly necessary for quality control of all germline datasets. For example, the gnomAD dataset comprises both exomes and genomes, and a substantial portion of the code is written to handle technical differences between those call sets, as well as to perform relevant joint analyses (such as inferring cryptically related individuals across exomes and genomes). These steps may not be relevant for all call sets.

Because of that, we have forked and modified the [gnomad_qc project](https://github.com/vladsaveliev/gnomad_qc), as well as the dependency project [gnomad_methods](https://github.com/vladsaveliev/gnomad_methods) to leave out parts referencing private Broad resources, drop the exome sample processing parts, and also parametrise some things like the output bucket and Slack notifications.
