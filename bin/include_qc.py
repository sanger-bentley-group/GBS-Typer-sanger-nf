#!/usr/bin/env python3
import argparse
import sys
import pandas as pd

from collections import defaultdict
from lib.get_headers import read_header_json
from lib.file_utils import FileUtils


def get_content(file, id, origin_names, new_names):

    df = pd.read_csv(file, sep="\t")
    origin_names.insert(0, id)
    df = df[origin_names]

    new_names.insert(0, id)
    df.columns = new_names

    return df


def create_report(qc_summary, qc_complete, in_silico, headers):

    qc_summary_df = get_content(qc_summary, headers["id"][0], list(headers["qc_summary"].keys()), list(headers["qc_summary"].values()))
    qc_complete_df = get_content(qc_complete, headers["id"][0], list(headers["qc_complete"].keys()), list(headers["qc_complete"].values()))
    in_silico_df = get_content(in_silico, headers["id"][0], list(headers["combine_all"].values()), list(headers["combine_all"].values()))

    qc_df = pd.merge(qc_summary_df, qc_complete_df, on = headers["id"], how = "outer")
    df = pd.merge(qc_df, in_silico_df, on = headers["id"], how = "outer")

    column_names = headers["id"]
    column_names.extend(list(headers["qc_summary"].values()))
    column_names.extend(list(headers["qc_complete"].values()))
    column_names.extend(list(headers["combine_all"].values()))

    model_df = pd.DataFrame(columns=column_names, index = [0])
    final_df = model_df.combine_first(df)

    return final_df


def get_arguments():

    parser = argparse.ArgumentParser(description='Combine QC reports and in silico report.')
    parser.add_argument('-s', '--qc_summary', dest='qc_summary', required=True, help='QC summary report.')
    parser.add_argument('-c', '--qc_complete', dest='qc_complete', required=True, help='QC complete report.')
    parser.add_argument('-j', '--headers', dest='header', required=True, help='Header file.')
    parser.add_argument('-i', '--in_silico', dest='in_silico', required=True, help='In silico report.')
    parser.add_argument('-o', '--output', dest='output', required=True, type=str, help='Output file.')

    return parser


def main(args):

    header_dict = read_header_json(args.header)

    df = create_report(args.qc_summary, args.qc_complete, args.in_silico, header_dict)

    FileUtils.write_pandas_output(df, args.output)


if __name__ == '__main__':
    parser = get_arguments()
    args = parser.parse_args()
    sys.exit(main(args))
