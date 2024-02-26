import polars as pl
from pathlib import Path
import pathlib
import sys

# Helper function to save individual maf files
def save_maf(
        incoming_maf,
        out_path = str
):
    # Account for the case when empty maf is passed
    try:
        name = incoming_maf['Tumor_Sample_Barcode'].unique().to_list()[0]
    except IndexError:
        return None
    path = f"{out_path}/{name}.maf"

    # Convert to a regular DataFrame and write to a TSV file
    incoming_maf.write_csv(path, separator='\t')

# Split multi-sample maf file to individual files
# at user-specified location
def prepare_sbs_mafs(
        out_path = str,
        subset_to_panel = False,
        panel = None,
        **mafs
):
    print('Preprocessing incoming maf file ...')

    if 'file_path' in mafs:
        maf_data = pl.read_csv(
            source = mafs['file_path'],
            has_header = True,
            separator = '\t',
            comment_prefix = '#',
            dtypes = {
                "Chromosome": pl.String,
                "Tumor_Sample_Barcode": pl.String
            }
        )

    elif 'maf_data' in mafs:
        maf_data = mafs['maf_data']

    else:
        print('Neither file_path nor maf_data is provided!')
        sys.exit('Please provide maf data as path to file or data frame')

    # Group by 'Tumor_Sample_Barcode' and split into polars DataFrames
    maf_data_grouped = maf_data.group_by('Tumor_Sample_Barcode')

    ind_mafs = [group for _, group in maf_data_grouped]

    Path(out_path).mkdir(parents = True, exist_ok = True)

    for maf in ind_mafs:
        if subset_to_panel:
            maf = cool_overlaps(maf, panel)
        save_maf(maf, out_path = out_path)

    return(ind_mafs)

# Estimate exposure based on maf files
from SigProfilerAssignment import Analyzer as Analyze

def run_sigprofiler(
        incoming_data = str,
        out_path = str,
        genome_build = "GRCh37"
):
    print('Running SigProfiler ...')
    Analyze.cosmic_fit(
        samples = incoming_data,
        output = out_path,
        input_type = "vcf",
        context_type = "96",
        genome_build = genome_build,
        cosmic_version = 3.4)
    return()

# Normalize signature exposure to be relative/sample
import polars.selectors as cs
def scale_sbs_exposure(
        file_path = str
):
    print('Scaling SBS exposure per sample ...')
    activities = pl.read_csv(
        source = file_path,
        has_header = True,
        separator = "\t",
        dtypes = {
            'Samples' : pl.String
        }
    )

    # Calculate sum of numeric columns horizontally
    # Store the sum in temp column `row_sum`
    activities = activities.with_columns(
        activities.select(
            cs.integer()
        ).sum(
            axis = 1
        ).alias(
            'row_sum'
        )
    )

    # List of column names
    numeric_cols = activities.select(cs.integer()).columns

    # Scale numeric values horizontally
    activities[numeric_cols] = activities[numeric_cols] / activities['row_sum']

    # Drop the temp column
    activities = activities.select(pl.col('*').exclude('row_sum'))

    return(activities)

def select_represented_sbs(
        input_data,
        threshold = 0.1
):
    # Only select columns with SBS data
    filtered = (
        input_data
        .select([col for col in input_data.columns if 'SBS' in col])
    )

    # Which SBS columns are not all 0s?
    sum_non_zero = filtered.select(pl.all().map(non_zero_percentage))

    # Only keep columns where there is at least X samples with non-0 SBS
    represented = sum_non_zero > (len(filtered) * threshold)

    columns = [col for col in represented.columns if represented[col][0] == True]

    represented = filtered[columns]

    # Add back the sample ids
    selected = represented.insert_at_idx(0, input_data['Samples'])

    return(represented)

import shutil
from .helpers import *
def estimate_sbs_exposure(
        out_path = str,
        genome_build = "GRCh37",
        clear_temp_outputs = True,
        subset_to_panel = False,
        panel = None,
        **mafs
):
    preprocess = prepare_sbs_mafs(
        out_path = out_path,
        subset_to_panel = subset_to_panel,
        panel = panel,
        **mafs
    )
    estimate = run_sigprofiler(
        incoming_data = out_path,
        out_path = out_path,
        genome_build = genome_build
    )
    input_file = sorted(pathlib.Path(out_path).glob('**/*_Activities.txt'))
    activities = scale_sbs_exposure(
        file_path = str(input_file[0])
    )

    if clear_temp_outputs:
        shutil.rmtree(out_path, ignore_errors = True)

    return(activities)
