# Helper functions
def non_zero_percentage(col):
    return col.filter(col > 0).count()

import polars as pl
def cool_overlaps(
        df1,
        df2,
        col_name1 = 'Chromosome',
        start_name1 = 'Start_Position',
        end_name1 = 'End_Position',
        col_name2 = 'chrom',
        start_name2 = 'start',
        end_name2 = 'end'
):

    # Ensure the kernel is not killed for large mafs
    # Use arbitrary cutoff of 50k variants
    if df1.shape[0] > 5e4:
        chromosomes = df1[col_name1].unique().to_list()
        chromosomes.sort()
        result_df = pl.DataFrame()
        for chromosome in chromosomes:
            small_maf = df1.filter(pl.col(col_name1) == chromosome)
            small_panel = df2.filter(pl.col(col_name2) == chromosome)
            overlap = small_maf.join(
                small_panel,
                left_on = col_name1,
                right_on = col_name2
            )

            overlap = overlap.filter(
                (overlap[start_name1] <= overlap[end_name2]) &
                (overlap[end_name1] >= overlap[start_name2])
            )

            overlap = overlap[small_maf.columns]
            result_df = pl.concat([result_df, overlap])

            return(result_df)

    result_df = df1.join(
        df2,
        left_on = col_name1,
        right_on = col_name2
    )

    overlap = result_df.filter(
        (result_df[start_name1] <= result_df[end_name2]) &
        (result_df[end_name1] >= result_df[start_name2])
    )

    overlap = overlap[df1.columns]

    return(overlap)

sbs_colors_list = [["SBS1", "#acf2d0"],
            ["SBS5", "#63d69e"],
            ["SBS2", "#f8b6b3"],
            ["SBS13", "#f17fb2"],
            ["SBS3", "#c4abc4"],
            ["SBS4", "#bcf2f5"],
            ["SBS7a", "#b5d7f5"],
            ["SBS7b", "#9ecef7"],
            ["SBS7c", "#84bdf0"],
            ["SBS7d", "#6cb2f0"],
            ["SBS8", "#dfc4f5"],
            ["SBS9", "#ebf5bc"],
            ["SBS10a", "#f2aeae"],
            ["SBS10b", "#f08080"],
            ["SBS17a", "#d9f7b0"],
            ["SBS17b", "#8cc63f"],
            ["SBS40", "#c4c4f5"],
            ["SBS6", "#faf1dc"],
            ["SBS14", "#faecca"],
            ["SBS15", "#fcebc2"],
            ["SBS20", "#fae4af"],
            ["SBS21", "#fae1a5"],
            ["SBS26", "#fcde97"],
            ["SBS44", "#fad682"],
            ["SBS27", "#C8C8C8"],
            ["SBS43", "#C0C0C0"],
            ["SBS45", "#BEBEBE"],
            ["SBS46", "#B8B8B8"],
            ["SBS47", "#B0B0B0"],
            ["SBS48", "#A9A9A9"],
            ["SBS49", "#A8A8A8"],
            ["SBS50", "#A0A0A0"],
            ["SBS51", "#989898"],
            ["SBS52", "#909090"],
            ["SBS53", "#888888"],
            ["SBS54", "#808080"],
            ["SBS55", "#787878"],
            ["SBS56", "#707070"],
            ["SBS57", "#696969"],
            ["SBS58", "#686868"],
            ["SBS59", "#606060"],
            ["SBS60", "#585858"],
            ["SBS37", "#EE7B73"],
            ["SBS84", "#486ad9"],
            ["SBS85", "#cd77d1"],
            ["SBS87", "#98ba93"],
            ["SBS19", "#bdbc45"],
            ["SBS40c", "#f5bd93"],
            ["SBS30", "#5f4ea1"],
            ["SBS93", "#709045"],
            ["SBS36", "#a3635b"],
            ["SBS32", "#376436"],
            ["SBS10c", "#ed6d6d"],
            ["SBS31", "#d56e37"],
            ["SBS16", "#826aaf"],
            ["SBS92", "#97ab74"],
            ["SBS10d", "#cc80af"],
            ["SBS98", "#d15a92"],
            ["SBS88", "#3f54a0"],
            ["SBS39", "#ad5a9c"],
            ["SBS34", "#b73181"],
            ["SBS35", "#c92487"],
            ["SBS96", "#85cbdb"],
            ["SBS95", "#80a5ad"],
            ["SBS33", "#5986c0"],
            ["SBS25", "#163816"],
            ["SBS99", "#603c90"],
            ["SBS23", "#529f52"],
            ["SBS40a", "#f5ddcb"],
            ["SBS28", "#da4090"],
            ["SBS12", "#fab9d6"],
            ["SBS29", "#de5634"],
            ["SBS22b", "#c73c35"],
            ["SBS97", "#e074a6"],
            ["SBS91", "#332d70"],
            ["SBS90", "#474088"],
            ["SBS86", "#492e78"],
            ["SBS41", "#9c9cf0"],
            ["SBS24", "#58bdcd"],
            ["SBS22a", "#85594f"],
            ["SBS18", "#bfbe5c"],
            ["SBS94", "#465c29"],
            ["SBS89", "#ce7693"],
            ["SBS42", "#478a4b"],
            ["SBS11", "#e6863c"],
            ["SBS38", "#ab679d"],
            ["SBS40b", "#eaa76c"]
]

pairing_colors_list = [
    ["matched", "#F0B67F"],
    ["unmatched", "#D6D1B1"],
    ["NA", "#ACADAF"]
]

ffpe_colors_list = [
    ["ff", "#009FFD"],
    ["FF", "#009FFD"],
    ["frozen", "#009FFD"],
    ["ffpe", "#95B2B8"],
    ["FFPE", "#95B2B8"],
    ["NA", "#ACADAF"]
]

pathology_colors_list = [
    ["BL", "#926CAD"],
    ["FL", "#EA8368"],
    ["DLBCL", "#479450"],
    ["COMFL", "#8BBC98"],
    ["PBL", "#E058C0"],
    ["B-ALL", "#C1C64B"],
    ["MZL", "#065A7F"],
    ["DLBCL-BL-like", "#34C7F4"],
    ["HGBL", "#294936"],
    ["NA", "#ACADAF"]
]

lymphgen_colors_list = [
    ["EZB-MYC", "#52000F"],
    ["EZB", "#721F0F"],
    ["EZB-COMP", "#C7371A"],
    ["ST2", "#C41230"],
    ["ST2-COMP", "#EC3251"],
    ["MCD", "#3B5FAC"],
    ["MCD-COMP", "#6787CB"],
    ["BN2", "#7F3293"],
    ["BN2-COMP", "#A949C1"],
    ["N1", "#55B55E"],
    ["N1-COMP", "#7FC787"],
    ["A53", "#5b6d8a"],
    ["Other", "#ACADAF"],
    ["COMPOSITE", "#ACADAF"],
    ["NA", "#ACADAF"]
]

colors = {
    'sbs': {item[0]: item[1] for item in sbs_colors_list},
    'pairing_status': {item[0]: item[1] for item in pairing_colors_list},
    'ffpe': {item[0]: item[1] for item in ffpe_colors_list},
    'pathology': {item[0]: item[1] for item in pathology_colors_list},
    'lymphgen': {item[0]: item[1] for item in lymphgen_colors_list}
}

maf_header = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
    "Chromosome", "Start_Position", "End_Position", "Strand",
    "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS",
    "dbSNP_Val_Status", "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1",
    "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
    "Match_Norm_Validation_Allele2", "Verification_Status",
    "Validation_Status", "Mutation_Status", "Sequencing_Phase",
    "Sequence_Source", "Validation_Method", "Score", "BAM_File",
    "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID",
    "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID",
    "Exon_Number", "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count"
]