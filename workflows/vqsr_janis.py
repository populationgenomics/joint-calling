"""
This is a conversion from the ukbiobank WDL workflow to Janis
"""
from datetime import datetime
from typing import List, Optional, Dict, Any

from janis_core import *
from janis_core.types.common_data_types import (
    String,
    File,
    Int,
    Array,
    Float,
    Boolean,
)

Tabixbgzippedfile_Dev = CommandToolBuilder(
    tool="TabixBGzippedFile",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="zipped_vcf", input_type=File(), doc=InputDocumentation(doc=None)
        ),
        ToolInput(
            tag="zipped_vcf_path",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="localized_tabix_path",
            input_type=String(),
            default=AddOperator(
                InputSelector(input_to_select="zipped_vcf", type_hint=File()),
                ".tbi",
            ),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="bucket_tabix_path",
            input_type=String(),
            default=AddOperator(
                InputSelector(input_to_select="zipped_vcf_path", type_hint=File()),
                ".tbi",
            ),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="bucket_tabix_output",
            output_type=String(),
            selector=InputSelector(
                input_to_select="bucket_tabix_path", type_hint=File()
            ),
            doc=OutputDocumentation(doc=None),
        )
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    gatk IndexFeatureFile -F {JANIS_WDL_TOKEN_1}\n    gsutil cp {JANIS_WDL_TOKEN_1}'.tbi' {JANIS_WDL_TOKEN_2}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="zipped_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="bucket_tabix_path", type_hint=File()
            ),
        )
    },
)
Splitintervallist_Dev = CommandToolBuilder(
    tool="SplitIntervalList",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="interval_list",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="scatter_count",
            input_type=Int(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="ref_fasta", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="ref_fasta_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="ref_dict", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="sample_names_unique_done",
            input_type=Boolean(),
            doc=InputDocumentation(doc=None, quality=InputQualityType.configuration),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="output_intervals",
            output_type=Array(t=File()),
            selector=WildcardSelector(wildcard="scatterDir/*"),
            doc=OutputDocumentation(doc=None),
        )
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    gatk --java-options -Xms3g SplitIntervals \\n      -L {JANIS_WDL_TOKEN_1} -O  scatterDir -scatter {JANIS_WDL_TOKEN_2} -R {JANIS_WDL_TOKEN_3} \\n      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW\n   ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="interval_list", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="scatter_count", type_hint=File()
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="ref_fasta", type_hint=File()
            ),
        )
    },
)
Gnarlygenotyperonvcf_Dev = CommandToolBuilder(
    tool="GnarlyGenotyperOnVcf",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="combined_gvcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="combined_gvcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="interval", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="output_vcf_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="ref_fasta", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="ref_fasta_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="ref_dict", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="dbsnp_vcf", input_type=String(), doc=InputDocumentation(doc=None)
        ),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="output_vcf",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="output_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="output_vcf_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.tbi",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="output_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -e\n\n    gatk --java-options -Xms8g \\n      GnarlyGenotyper \\n      -R {JANIS_WDL_TOKEN_1} \\n      -O {JANIS_WDL_TOKEN_2} \\n      -D {JANIS_WDL_TOKEN_3} \\n      --only-output-calls-starting-in-intervals \\n      --keep-all-sites \\n      -V {JANIS_WDL_TOKEN_4} \\n      -L {JANIS_WDL_TOKEN_5}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="ref_fasta", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="output_vcf_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="dbsnp_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="combined_gvcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="interval", type_hint=File()
            ),
        )
    },
)
Hardfilterandmakesitesonlyvcf_Dev = CommandToolBuilder(
    tool="HardFilterAndMakeSitesOnlyVcf",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(tag="vcf", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(tag="vcf_index", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="excess_het_threshold",
            input_type=Float(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="variant_filtered_vcf_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_vcf_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="variant_filtered_vcf",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="variant_filtered_vcf_filename",
                    type_hint=File(),
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="variant_filtered_vcf_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.tbi",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="variant_filtered_vcf_filename",
                    type_hint=File(),
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="sites_only_vcf",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="sites_only_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="sites_only_vcf_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.tbi",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="sites_only_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    gatk --java-options -Xms3g \\n      VariantFiltration \\n      --filter-expression 'ExcessHet > {JANIS_WDL_TOKEN_1}' \\n      --filter-name ExcessHet \\n      -O {JANIS_WDL_TOKEN_2} \\n      -V {JANIS_WDL_TOKEN_3}\n\n    gatk --java-options -Xms3g \\n      MakeSitesOnlyVcf \\n      -I {JANIS_WDL_TOKEN_2} \\n      -O {JANIS_WDL_TOKEN_4}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="excess_het_threshold", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="variant_filtered_vcf_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_3=InputSelector(input_to_select="vcf", type_hint=File()),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="sites_only_vcf_filename", type_hint=File()
            ),
        )
    },
)
Gathervcfs_Dev = CommandToolBuilder(
    tool="GatherVcfs",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="input_vcfs",
            input_type=Array(t=File()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="output_vcf_name",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="output_vcf",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="output_vcf_name", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="output_vcf_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.tbi",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="output_vcf_name", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.\n    # This argument disables expensive checks that the file headers contain the same set of\n    # genotyped samples and that files are in order by position of first record.\n    gatk --java-options -Xms6g \\n      GatherVcfsCloud \\n      --ignore-safety-checks \\n      --gather-type BLOCK \\n      --input {JANIS_WDL_TOKEN_1} \\n      --output {JANIS_WDL_TOKEN_2}\n\n    tabix {JANIS_WDL_TOKEN_2}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="input_vcfs", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="output_vcf_name", type_hint=File()
            ),
        )
    },
)
Indelsvariantrecalibrator_Dev = CommandToolBuilder(
    tool="IndelsVariantRecalibrator",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="recalibration_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="tranches_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="model_report",
            input_type=File(optional=True),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_tranche_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_annotation_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="mills_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="axiomPoly_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="mills_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="axiomPoly_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="use_allele_specific_annotations",
            input_type=Boolean(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="max_gaussians",
            input_type=Int(),
            default=4,
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
    ],
    outputs=[
        ToolOutput(
            tag="recalibration",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibration_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="recalibration_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.idx",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibration_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="tranches",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="tranches_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="ubuntu:latest",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    MODEL_REPORT={JANIS_WDL_TOKEN_1}\n\n    gatk --java-options -Xms100g \\n      VariantRecalibrator \\n      -V {JANIS_WDL_TOKEN_2} \\n      -O {JANIS_WDL_TOKEN_3} \\n      --tranches-file {JANIS_WDL_TOKEN_4} \\n      --trust-all-polymorphic \\n      -tranche {JANIS_WDL_TOKEN_5} \\n      -an {JANIS_WDL_TOKEN_6} \\n      -mode INDEL \\n      {JANIS_WDL_TOKEN_7} \\n      {JANIS_WDL_TOKEN_8} \\n      --max-gaussians {JANIS_WDL_TOKEN_9} \\n      -resource:mills,known=false,training=true,truth=true,prior=12 {JANIS_WDL_TOKEN_10} \\n      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {JANIS_WDL_TOKEN_11} \\n      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {JANIS_WDL_TOKEN_12}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="model_report", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="sites_only_variant_filtered_vcf",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="recalibration_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="tranches_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="recalibration_tranche_values", type_hint=File()
            ),
            JANIS_WDL_TOKEN_6=InputSelector(
                input_to_select="recalibration_annotation_values",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_7=InputSelector(
                input_to_select="use_allele_specific_annotations",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_8=InputSelector(
                input_to_select="model_report_arg", type_hint=File()
            ),
            JANIS_WDL_TOKEN_9=InputSelector(
                input_to_select="max_gaussians", type_hint=File()
            ),
            JANIS_WDL_TOKEN_10=InputSelector(
                input_to_select="mills_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_11=InputSelector(
                input_to_select="axiomPoly_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_12=InputSelector(
                input_to_select="dbsnp_resource_vcf", type_hint=File()
            ),
        )
    },
)
Snpsvariantrecalibratorcreatemodel_Dev = CommandToolBuilder(
    tool="SNPsVariantRecalibratorCreateModel",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="recalibration_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="tranches_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="downsampleFactor",
            input_type=Int(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="model_report_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_tranche_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_annotation_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="hapmap_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="omni_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="one_thousand_genomes_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="hapmap_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="omni_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="one_thousand_genomes_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="max_gaussians",
            input_type=Int(),
            default=6,
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="java_mem",
            input_type=Int(),
            default=100,
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="use_allele_specific_annotations",
            input_type=Boolean(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.4.1",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="model_report",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="model_report_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        )
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.4.1",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    gatk --java-options -Xms{JANIS_WDL_TOKEN_1}g \\n      VariantRecalibrator \\n      -V {JANIS_WDL_TOKEN_2} \\n      -O {JANIS_WDL_TOKEN_3} \\n      --tranches-file {JANIS_WDL_TOKEN_4} \\n      --trust-all-polymorphic \\n      -tranche {JANIS_WDL_TOKEN_5} \\n      -an {JANIS_WDL_TOKEN_6} \\n      -mode SNP \\n      {JANIS_WDL_TOKEN_7} \\n      --sample-every-Nth-variant {JANIS_WDL_TOKEN_8} \\n      --output-model {JANIS_WDL_TOKEN_9} \\n      --max-gaussians {JANIS_WDL_TOKEN_10} \\n      -resource:hapmap,known=false,training=true,truth=true,prior=15 {JANIS_WDL_TOKEN_11} \\n      -resource:omni,known=false,training=true,truth=true,prior=12 {JANIS_WDL_TOKEN_12} \\n      -resource:1000G,known=false,training=true,truth=false,prior=10 {JANIS_WDL_TOKEN_13} \\n      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {JANIS_WDL_TOKEN_14}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="java_mem", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="sites_only_variant_filtered_vcf",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="recalibration_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="tranches_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="recalibration_tranche_values", type_hint=File()
            ),
            JANIS_WDL_TOKEN_6=InputSelector(
                input_to_select="recalibration_annotation_values",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_7=InputSelector(
                input_to_select="use_allele_specific_annotations",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_8=InputSelector(
                input_to_select="downsampleFactor", type_hint=File()
            ),
            JANIS_WDL_TOKEN_9=InputSelector(
                input_to_select="model_report_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_10=InputSelector(
                input_to_select="max_gaussians", type_hint=File()
            ),
            JANIS_WDL_TOKEN_11=InputSelector(
                input_to_select="hapmap_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_12=InputSelector(
                input_to_select="omni_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_13=InputSelector(
                input_to_select="one_thousand_genomes_resource_vcf",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_14=InputSelector(
                input_to_select="dbsnp_resource_vcf", type_hint=File()
            ),
        )
    },
)
Snpsvariantrecalibrator_Dev = CommandToolBuilder(
    tool="SNPsVariantRecalibrator",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="recalibration_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="tranches_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="model_report",
            input_type=File(optional=True),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_tranche_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="recalibration_annotation_values",
            input_type=Array(t=String()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="sites_only_variant_filtered_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="hapmap_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="omni_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="one_thousand_genomes_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="hapmap_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="omni_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="one_thousand_genomes_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="dbsnp_resource_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="max_gaussians",
            input_type=Int(),
            default=6,
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="machine_mem_gb",
            input_type=Int(optional=True),
            doc=InputDocumentation(doc=None, quality=InputQualityType.configuration),
        ),
        ToolInput(
            tag="use_allele_specific_annotations",
            input_type=Boolean(),
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="recalibration",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibration_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="recalibration_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.idx",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibration_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="tranches",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="tranches_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    MODEL_REPORT={JANIS_WDL_TOKEN_1}\n\n    gatk --java-options -Xms{JANIS_WDL_TOKEN_2}g \\n      VariantRecalibrator \\n      -V {JANIS_WDL_TOKEN_3} \\n      -O {JANIS_WDL_TOKEN_4} \\n      --tranches-file {JANIS_WDL_TOKEN_5} \\n      --trust-all-polymorphic \\n      -tranche {JANIS_WDL_TOKEN_6} \\n      -an {JANIS_WDL_TOKEN_7} \\n      -mode SNP \\n      {JANIS_WDL_TOKEN_8} \\n       {JANIS_WDL_TOKEN_9} \\n      --max-gaussians {JANIS_WDL_TOKEN_10} \\n      -resource:hapmap,known=false,training=true,truth=true,prior=15 {JANIS_WDL_TOKEN_11} \\n      -resource:omni,known=false,training=true,truth=true,prior=12 {JANIS_WDL_TOKEN_12} \\n      -resource:1000G,known=false,training=true,truth=false,prior=10 {JANIS_WDL_TOKEN_13} \\n      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {JANIS_WDL_TOKEN_14}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="model_report", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=SubtractOperator(
                InputSelector(input_to_select="machine_mem", type_hint=File()), 1
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="sites_only_variant_filtered_vcf",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="recalibration_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="tranches_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_6=InputSelector(
                input_to_select="recalibration_tranche_values", type_hint=File()
            ),
            JANIS_WDL_TOKEN_7=InputSelector(
                input_to_select="recalibration_annotation_values",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_8=InputSelector(
                input_to_select="use_allele_specific_annotations",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_9=InputSelector(
                input_to_select="model_report_arg", type_hint=File()
            ),
            JANIS_WDL_TOKEN_10=InputSelector(
                input_to_select="max_gaussians", type_hint=File()
            ),
            JANIS_WDL_TOKEN_11=InputSelector(
                input_to_select="hapmap_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_12=InputSelector(
                input_to_select="omni_resource_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_13=InputSelector(
                input_to_select="one_thousand_genomes_resource_vcf",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_14=InputSelector(
                input_to_select="dbsnp_resource_vcf", type_hint=File()
            ),
        )
    },
)
Gathertranches_Dev = CommandToolBuilder(
    tool="GatherTranches",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="tranches",
            input_type=Array(t=File()),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="output_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="out_tranches",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="output_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        )
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    tranches_fofn={JANIS_WDL_TOKEN_1}\n\n    # Jose says:\n    # Cromwell will fall over if we have it try to localize tens of thousands of files,\n    # so we manually localize files using gsutil.\n    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)\n    # PAPI doesn't do.\n\n    # This is here to deal with the JES bug where commands may be run twice\n    rm -rf tranches\n    mkdir tranches\n    RETRY_LIMIT=5\n\n    count=0\n    until cat $tranches_fofn | /usr/bin/gsutil -m cp -L cp.log -c -I tranches/; do\n        sleep 1\n        ((count++)) && ((count >= $RETRY_LIMIT)) && break\n    done\n    if [ '$count' -ge '$RETRY_LIMIT' ]; then\n        echo 'Could not copy all the tranches from the cloud' && exit 1\n    fi\n\n    cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print 'tranches/' $1}' > inputs.list\n\n    gatk --java-options -Xms6g \\n      GatherTranches \\n      --input inputs.list \\n      --output {JANIS_WDL_TOKEN_2}\n  ",
            JANIS_WDL_TOKEN_1="JANIS: write_lines([inputs.tranches])",
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="output_filename", type_hint=File()
            ),
        )
    },
)
Applyrecalibration_Dev = CommandToolBuilder(
    tool="ApplyRecalibration",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(
            tag="recalibrated_vcf_filename",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="input_vcf", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="input_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="indels_recalibration",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="indels_recalibration_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="indels_tranches",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="snps_recalibration",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="snps_recalibration_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="snps_tranches",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="indel_filter_level",
            input_type=Float(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="snp_filter_level",
            input_type=Float(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="use_allele_specific_annotations",
            input_type=Boolean(),
            doc=InputDocumentation(doc=None, quality=InputQualityType.configuration),
        ),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="recalibrated_vcf",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibrated_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="recalibrated_vcf_index",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.tbi",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="recalibrated_vcf_filename", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    gatk --java-options -Xms5g \\n      ApplyVQSR \\n      -O tmp.indel.recalibrated.vcf \\n      -V {JANIS_WDL_TOKEN_1} \\n      --recal-file {JANIS_WDL_TOKEN_2} \\n      --tranches-file {JANIS_WDL_TOKEN_3} \\n      --truth-sensitivity-filter-level {JANIS_WDL_TOKEN_4} \\n      --create-output-variant-index true \\n      -mode INDEL {JANIS_WDL_TOKEN_5} \\n\n\n    gatk --java-options -Xms5g \\n      ApplyVQSR \\n      -O {JANIS_WDL_TOKEN_6} \\n      -V tmp.indel.recalibrated.vcf \\n      --recal-file {JANIS_WDL_TOKEN_7} \\n      --tranches-file {JANIS_WDL_TOKEN_8} \\n      --truth-sensitivity-filter-level {JANIS_WDL_TOKEN_9} \\n      --create-output-variant-index true \\n      -mode SNP {JANIS_WDL_TOKEN_5} \\n\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="input_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="indels_recalibration", type_hint=File()
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="indels_tranches", type_hint=File()
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="indel_filter_level", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="use_allele_specific_annotations",
                type_hint=File(),
            ),
            JANIS_WDL_TOKEN_6=InputSelector(
                input_to_select="recalibrated_vcf_filename", type_hint=File()
            ),
            JANIS_WDL_TOKEN_7=InputSelector(
                input_to_select="snps_recalibration", type_hint=File()
            ),
            JANIS_WDL_TOKEN_8=InputSelector(
                input_to_select="snps_tranches", type_hint=File()
            ),
            JANIS_WDL_TOKEN_9=InputSelector(
                input_to_select="snp_filter_level", type_hint=File()
            ),
        )
    },
)
Collectvariantcallingmetrics_Dev = CommandToolBuilder(
    tool="CollectVariantCallingMetrics",
    base_command=["sh", "script.sh"],
    inputs=[
        ToolInput(tag="input_vcf", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="input_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="metrics_filename_prefix",
            input_type=String(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="dbsnp_vcf", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="dbsnp_vcf_index",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="interval_list",
            input_type=File(),
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(tag="ref_dict", input_type=File(), doc=InputDocumentation(doc=None)),
        ToolInput(tag="disk_size", input_type=Int(), doc=InputDocumentation(doc=None)),
        ToolInput(
            tag="gatk_docker",
            input_type=String(),
            default="us.gcr.io/broad-gatk/gatk:4.1.1.0",
            doc=InputDocumentation(doc=None),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="detail_metrics_file",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.variant_calling_detail_metrics",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="metrics_filename_prefix", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
        ToolOutput(
            tag="summary_metrics_file",
            output_type=File(),
            selector=StringFormatter(
                "{JANIS_WDL_TOKEN_1}.variant_calling_summary_metrics",
                JANIS_WDL_TOKEN_1=InputSelector(
                    input_to_select="metrics_filename_prefix", type_hint=File()
                ),
            ),
            doc=OutputDocumentation(doc=None),
        ),
    ],
    container="us.gcr.io/broad-gatk/gatk:4.1.1.0",
    version="DEV",
    files_to_create={
        "script.sh": StringFormatter(
            "\n    set -euo pipefail\n\n    gatk --java-options -Xms6g \\n      CollectVariantCallingMetrics \\n      --INPUT {JANIS_WDL_TOKEN_1} \\n      --DBSNP {JANIS_WDL_TOKEN_2} \\n      --SEQUENCE_DICTIONARY {JANIS_WDL_TOKEN_3} \\n      --OUTPUT {JANIS_WDL_TOKEN_4} \\n      --THREAD_COUNT 8 \\n      --TARGET_INTERVALS {JANIS_WDL_TOKEN_5}\n  ",
            JANIS_WDL_TOKEN_1=InputSelector(
                input_to_select="input_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_2=InputSelector(
                input_to_select="dbsnp_vcf", type_hint=File()
            ),
            JANIS_WDL_TOKEN_3=InputSelector(
                input_to_select="ref_dict", type_hint=File()
            ),
            JANIS_WDL_TOKEN_4=InputSelector(
                input_to_select="metrics_filename_prefix", type_hint=File()
            ),
            JANIS_WDL_TOKEN_5=InputSelector(
                input_to_select="interval_list", type_hint=File()
            ),
        )
    },
)

Variantcallingofthefuture = WorkflowBuilder(
    identifier="VariantCallingOFTHEFUTURE",
)

Variantcallingofthefuture.input(
    "unpadded_intervals_file",
    String(),
)

Variantcallingofthefuture.input(
    "combined_gvcf",
    File(),
)

Variantcallingofthefuture.input(
    "callset_name",
    String(),
)

Variantcallingofthefuture.input(
    "ref_fasta",
    File(),
)

Variantcallingofthefuture.input(
    "ref_fasta_index",
    File(),
)

Variantcallingofthefuture.input(
    "ref_dict",
    File(),
)

Variantcallingofthefuture.input(
    "dbsnp_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "dbsnp_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "small_disk",
    Int(),
)

Variantcallingofthefuture.input(
    "medium_disk",
    Int(),
)

Variantcallingofthefuture.input(
    "huge_disk",
    Int(),
)

Variantcallingofthefuture.input(
    "snp_recalibration_tranche_values",
    Array(t=String()),
)

Variantcallingofthefuture.input(
    "snp_recalibration_annotation_values",
    Array(t=String()),
)

Variantcallingofthefuture.input(
    "indel_recalibration_tranche_values",
    Array(t=String()),
)

Variantcallingofthefuture.input(
    "indel_recalibration_annotation_values",
    Array(t=String()),
)

Variantcallingofthefuture.input(
    "eval_interval_list",
    File(),
)

Variantcallingofthefuture.input(
    "hapmap_resource_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "hapmap_resource_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "omni_resource_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "omni_resource_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "one_thousand_genomes_resource_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "one_thousand_genomes_resource_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "mills_resource_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "mills_resource_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "axiomPoly_resource_vcf",
    File(),
)

Variantcallingofthefuture.input(
    "axiomPoly_resource_vcf_index",
    File(),
)

Variantcallingofthefuture.input(
    "dbsnp_resource_vcf",
    File(optional=True),
    default=Variantcallingofthefuture.dbsnp_vcf,
)

Variantcallingofthefuture.input(
    "dbsnp_resource_vcf_index",
    File(optional=True),
    default=Variantcallingofthefuture.dbsnp_vcf_index,
)

Variantcallingofthefuture.input(
    "excess_het_threshold",
    Float(optional=True),
    default=54.69,
)

Variantcallingofthefuture.input(
    "snp_filter_level",
    Float(),
)

Variantcallingofthefuture.input(
    "indel_filter_level",
    Float(),
)

Variantcallingofthefuture.input(
    "SNP_VQSR_downsampleFactor",
    Int(),
)

Variantcallingofthefuture.input(
    "indel_VQSR_downsampleFactor",
    Int(),
)

Variantcallingofthefuture.input(
    "use_allele_specific_annotations",
    Boolean(optional=True),
    default=True,
)

Variantcallingofthefuture.input(
    "vcf_count",
    Int(),
)

Variantcallingofthefuture.input(
    "do_tabix",
    Boolean(),
)

Variantcallingofthefuture.input(
    "unboundedScatterCount",
    Int(optional=True),
    default=Variantcallingofthefuture.vcf_count,
)

Variantcallingofthefuture.input(
    "scatterCount",
    Int(optional=True),
    default=If(
        GtOperator(Variantcallingofthefuture.unboundedScatterCount, 10),
        Variantcallingofthefuture.unboundedScatterCount,
        10,
    ),
)

Variantcallingofthefuture.input(
    "SplitIntervalList_sample_names_unique_done",
    Boolean(optional=True),
    default=True,
)

Variantcallingofthefuture.input(
    "unpadded_intervals",
    Array(t=File(), optional=True),
)

Variantcallingofthefuture.input(
    "SNPsVariantRecalibratorScattered_machine_mem_gb",
    Int(optional=True),
    default=60,
)

Variantcallingofthefuture.input(
    "is_small_callset",
    Boolean(optional=True),
    default=False,
)

Variantcallingofthefuture.input(
    "ApplyRecalibration_use_allele_specific_annotations",
    Boolean(optional=True),
    default=True,
)


Variantcallingofthefuture.step(
    "TabixBGzippedFile",
    Tabixbgzippedfile_Dev(
        zipped_vcf=Variantcallingofthefuture.combined_gvcf,
        zipped_vcf_path=Variantcallingofthefuture.combined_gvcf,
    ),
    when=Variantcallingofthefuture.do_tabix,
)


Variantcallingofthefuture.step(
    "SplitIntervalList",
    Splitintervallist_Dev(
        interval_list=Variantcallingofthefuture.unpadded_intervals_file,
        scatter_count=Variantcallingofthefuture.scatterCount,
        ref_fasta=Variantcallingofthefuture.ref_fasta,
        ref_fasta_index=Variantcallingofthefuture.ref_fasta_index,
        ref_dict=Variantcallingofthefuture.ref_dict,
        disk_size=Variantcallingofthefuture.small_disk,
        sample_names_unique_done=Variantcallingofthefuture.SplitIntervalList_sample_names_unique_done,
    ),
)


Variantcallingofthefuture.step(
    "GnarlyGenotyperOnVcf",
    Gnarlygenotyperonvcf_Dev(
        combined_gvcf=Variantcallingofthefuture.combined_gvcf,
        combined_gvcf_index=FilterNullOperator(
            [
                Variantcallingofthefuture.TabixBGzippedFile.bucket_tabix_output,
                AddOperator(Variantcallingofthefuture.combined_gvcf, ".tbi"),
            ]
        ),
        interval=IndexOperator(
            Variantcallingofthefuture.unpadded_intervals, ForEachSelector()
        ),
        output_vcf_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, "."),
                ForEachSelector(),
            ),
            ".vcf.gz",
        ),
        ref_fasta=Variantcallingofthefuture.ref_fasta,
        ref_fasta_index=Variantcallingofthefuture.ref_fasta_index,
        ref_dict=Variantcallingofthefuture.ref_dict,
        dbsnp_vcf=Variantcallingofthefuture.dbsnp_vcf,
    ),
    foreach=RangeOperator(LengthOperator(Variantcallingofthefuture.unpadded_intervals)),
)


Variantcallingofthefuture.step(
    "HardFilterAndMakeSitesOnlyVcf",
    Hardfilterandmakesitesonlyvcf_Dev(
        vcf=Variantcallingofthefuture.GnarlyGenotyperOnVcf.output_vcf,
        vcf_index=Variantcallingofthefuture.GnarlyGenotyperOnVcf.output_vcf_index,
        excess_het_threshold=Variantcallingofthefuture.excess_het_threshold,
        variant_filtered_vcf_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, "."),
                ForEachSelector(),
            ),
            ".variant_filtered.vcf.gz",
        ),
        sites_only_vcf_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, "."),
                ForEachSelector(),
            ),
            ".sites_only.variant_filtered.vcf.gz",
        ),
        disk_size=Variantcallingofthefuture.medium_disk,
    ),
    foreach=RangeOperator(LengthOperator(Variantcallingofthefuture.unpadded_intervals)),
)


Variantcallingofthefuture.step(
    "SitesOnlyGatherVcf",
    Gathervcfs_Dev(
        input_vcfs=[
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.sites_only_vcf
        ],
        output_vcf_name=AddOperator(
            Variantcallingofthefuture.callset_name, ".sites_only.vcf.gz"
        ),
        disk_size=Variantcallingofthefuture.medium_disk,
    ),
)


Variantcallingofthefuture.step(
    "IndelsVariantRecalibrator",
    Indelsvariantrecalibrator_Dev(
        sites_only_variant_filtered_vcf=Variantcallingofthefuture.SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index=Variantcallingofthefuture.SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".indels.recal"
        ),
        tranches_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".indels.tranches"
        ),
        recalibration_tranche_values=Variantcallingofthefuture.indel_recalibration_tranche_values,
        recalibration_annotation_values=Variantcallingofthefuture.indel_recalibration_annotation_values,
        mills_resource_vcf=Variantcallingofthefuture.mills_resource_vcf,
        mills_resource_vcf_index=Variantcallingofthefuture.mills_resource_vcf_index,
        axiomPoly_resource_vcf=Variantcallingofthefuture.axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index=Variantcallingofthefuture.axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf=Variantcallingofthefuture.dbsnp_resource_vcf,
        dbsnp_resource_vcf_index=Variantcallingofthefuture.dbsnp_resource_vcf_index,
        use_allele_specific_annotations=Variantcallingofthefuture.use_allele_specific_annotations,
        disk_size=Variantcallingofthefuture.small_disk,
    ),
)


Variantcallingofthefuture.step(
    "SNPsVariantRecalibratorCreateModel",
    Snpsvariantrecalibratorcreatemodel_Dev(
        sites_only_variant_filtered_vcf=Variantcallingofthefuture.SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index=Variantcallingofthefuture.SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".snps.recal"
        ),
        tranches_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".snps.tranches"
        ),
        recalibration_tranche_values=Variantcallingofthefuture.snp_recalibration_tranche_values,
        recalibration_annotation_values=Variantcallingofthefuture.snp_recalibration_annotation_values,
        downsampleFactor=Variantcallingofthefuture.SNP_VQSR_downsampleFactor,
        model_report_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".snps.model.report"
        ),
        hapmap_resource_vcf=Variantcallingofthefuture.hapmap_resource_vcf,
        hapmap_resource_vcf_index=Variantcallingofthefuture.hapmap_resource_vcf_index,
        omni_resource_vcf=Variantcallingofthefuture.omni_resource_vcf,
        omni_resource_vcf_index=Variantcallingofthefuture.omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf=Variantcallingofthefuture.one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index=Variantcallingofthefuture.one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf=Variantcallingofthefuture.dbsnp_resource_vcf,
        dbsnp_resource_vcf_index=Variantcallingofthefuture.dbsnp_resource_vcf_index,
        disk_size=Variantcallingofthefuture.small_disk,
        use_allele_specific_annotations=Variantcallingofthefuture.use_allele_specific_annotations,
    ),
)


Variantcallingofthefuture.step(
    "SNPsVariantRecalibratorScattered",
    Snpsvariantrecalibrator_Dev(
        sites_only_variant_filtered_vcf=IndexOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
            ForEachSelector(),
        ),
        sites_only_variant_filtered_vcf_index=IndexOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
            ForEachSelector(),
        ),
        recalibration_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, ".snps."),
                ForEachSelector(),
            ),
            ".recal",
        ),
        tranches_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, ".snps."),
                ForEachSelector(),
            ),
            ".tranches",
        ),
        recalibration_tranche_values=Variantcallingofthefuture.snp_recalibration_tranche_values,
        recalibration_annotation_values=Variantcallingofthefuture.snp_recalibration_annotation_values,
        model_report=Variantcallingofthefuture.SNPsVariantRecalibratorCreateModel.model_report,
        hapmap_resource_vcf=Variantcallingofthefuture.hapmap_resource_vcf,
        hapmap_resource_vcf_index=Variantcallingofthefuture.hapmap_resource_vcf_index,
        omni_resource_vcf=Variantcallingofthefuture.omni_resource_vcf,
        omni_resource_vcf_index=Variantcallingofthefuture.omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf=Variantcallingofthefuture.one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index=Variantcallingofthefuture.one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf=Variantcallingofthefuture.dbsnp_resource_vcf,
        dbsnp_resource_vcf_index=Variantcallingofthefuture.dbsnp_resource_vcf_index,
        disk_size=Variantcallingofthefuture.small_disk,
        machine_mem_gb=Variantcallingofthefuture.SNPsVariantRecalibratorScattered_machine_mem_gb,
        use_allele_specific_annotations=Variantcallingofthefuture.use_allele_specific_annotations,
    ),
    foreach=RangeOperator(
        LengthOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.sites_only_vcf
        )
    ),
)


Variantcallingofthefuture.step(
    "SNPGatherTranches",
    Gathertranches_Dev(
        tranches=[Variantcallingofthefuture.SNPsVariantRecalibratorScattered.tranches],
        output_filename=AddOperator(
            Variantcallingofthefuture.callset_name, ".snps.gathered.tranches"
        ),
        disk_size=Variantcallingofthefuture.small_disk,
    ),
)


Variantcallingofthefuture.step(
    "ApplyRecalibration",
    Applyrecalibration_Dev(
        recalibrated_vcf_filename=AddOperator(
            AddOperator(
                AddOperator(Variantcallingofthefuture.callset_name, ".filtered."),
                ForEachSelector(),
            ),
            ".vcf.gz",
        ),
        input_vcf=IndexOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
            ForEachSelector(),
        ),
        input_vcf_index=IndexOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
            ForEachSelector(),
        ),
        indels_recalibration=Variantcallingofthefuture.IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index=Variantcallingofthefuture.IndelsVariantRecalibrator.recalibration_index,
        indels_tranches=Variantcallingofthefuture.IndelsVariantRecalibrator.tranches,
        snps_recalibration=IndexOperator(
            Variantcallingofthefuture.SNPsVariantRecalibratorScattered.recalibration,
            ForEachSelector(),
        ),
        snps_recalibration_index=IndexOperator(
            Variantcallingofthefuture.SNPsVariantRecalibratorScattered.recalibration_index,
            ForEachSelector(),
        ),
        snps_tranches=Variantcallingofthefuture.SNPGatherTranches.out_tranches,
        indel_filter_level=Variantcallingofthefuture.indel_filter_level,
        snp_filter_level=Variantcallingofthefuture.snp_filter_level,
        disk_size=Variantcallingofthefuture.medium_disk,
        use_allele_specific_annotations=Variantcallingofthefuture.ApplyRecalibration_use_allele_specific_annotations,
    ),
    foreach=RangeOperator(
        LengthOperator(
            Variantcallingofthefuture.HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf
        )
    ),
)


Variantcallingofthefuture.step(
    "CollectMetricsSharded",
    Collectvariantcallingmetrics_Dev(
        input_vcf=Variantcallingofthefuture.ApplyRecalibration.recalibrated_vcf,
        input_vcf_index=Variantcallingofthefuture.ApplyRecalibration.recalibrated_vcf_index,
        metrics_filename_prefix=AddOperator(
            AddOperator(Variantcallingofthefuture.callset_name, "."),
            ForEachSelector(),
        ),
        dbsnp_vcf=Variantcallingofthefuture.dbsnp_vcf,
        dbsnp_vcf_index=Variantcallingofthefuture.dbsnp_vcf_index,
        interval_list=Variantcallingofthefuture.eval_interval_list,
        ref_dict=Variantcallingofthefuture.ref_dict,
        disk_size=Variantcallingofthefuture.small_disk,
    ),
    when=NotOperator(Variantcallingofthefuture.is_small_callset),
)


Variantcallingofthefuture.step(
    "FinalGatherVcf",
    Gathervcfs_Dev(
        input_vcfs=[Variantcallingofthefuture.ApplyRecalibration.recalibrated_vcf],
        output_vcf_name=AddOperator(Variantcallingofthefuture.callset_name, ".vcf.gz"),
        disk_size=Variantcallingofthefuture.huge_disk,
    ),
    when=Variantcallingofthefuture.is_small_callset,
)


Variantcallingofthefuture.step(
    "CollectMetricsOnFullVcf",
    Collectvariantcallingmetrics_Dev(
        input_vcf=Variantcallingofthefuture.FinalGatherVcf.output_vcf,
        input_vcf_index=Variantcallingofthefuture.FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix=Variantcallingofthefuture.callset_name,
        dbsnp_vcf=Variantcallingofthefuture.dbsnp_vcf,
        dbsnp_vcf_index=Variantcallingofthefuture.dbsnp_vcf_index,
        interval_list=Variantcallingofthefuture.eval_interval_list,
        ref_dict=Variantcallingofthefuture.ref_dict,
        disk_size=Variantcallingofthefuture.huge_disk,
    ),
    when=Variantcallingofthefuture.is_small_callset,
)


if __name__ == "__main__":
    # or "cwl"
    from janis_core.translations.hailbatch import HailBatchTranslator

    tool = Variantcallingofthefuture()  # .translate("wdl")
    s = HailBatchTranslator.translate_workflow(tool)
    print(s)
