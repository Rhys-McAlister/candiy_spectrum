import os
import requests
import argparse
import logging
import pandas as pd
import time

from model.utils import set_logger

nist_url = "https://webbook.nist.gov/cgi/cbook.cgi"


def make_request_with_retry(url, params, max_retries=3, timeout=30):
    retry_count = 0
    while retry_count < max_retries:
        try:
            response = requests.get(url, params=params, timeout=timeout)
            return response
        except requests.exceptions.Timeout as e:
            logging.warning(
                f"Timeout error: {e}. Retrying {retry_count+1}/{max_retries}."
            )
            time.sleep(5)  # wait for 5 seconds before retrying
        except requests.exceptions.RequestException as e:
            logging.error(f"Request error: {e}. Skipping to next CAS number.")
            return None
        retry_count += 1
    logging.error(
        "Failed to get response after multiple retries. Skipping to next CAS number."
    )
    return None


def scrap_data(cas_ls, params, data_dir):
    spectra_path = os.path.join(data_dir, params["Type"].lower(), "")
    if not os.path.exists(spectra_path):
        os.makedirs(spectra_path)

    num_created = 0
    for cas_id in cas_ls:
        filename = os.path.join(spectra_path, cas_id + ".jdx")
        if os.path.exists(filename):
            logging.info(f"Skipping {filename}, already exists.")
            continue

        params["JCAMP"] = "C" + cas_id
        response = make_request_with_retry(nist_url, params)
        if response is None:
            continue

        if response.text == "##TITLE=Spectrum not found.\n##END=\n":
            continue

        num_created += 1
        logging.info(
            f'Creating {params["Type"].lower()} spectra for id: {cas_id}. Total spectra created {num_created}'
        )
        with open(filename, "wb") as data:
            data.write(response.content)


# def scrap_inchi(cas_ls, params, data_dir):
#     inchi_path = os.path.join(data_dir, "inchi.txt")
#     num_created = 0

#     existing_cas_ids = set()
#     if os.path.exists(inchi_path):
#         with open(inchi_path, "r") as file:
#             for line in file:
#                 cas_id = line.strip().split("\t")
#                 existing_cas_ids.add(cas_id)

#     with open(inchi_path, "a") as file:
#         if not existing_cas_ids:
#             file.write("cas_id\tinchi\n")

#         for cas_id in cas_ls:
#             if cas_id in existing_cas_ids:
#                 logging.info(f"Skipping InChi key for id: {cas_id}, already exists.")
#                 continue

#             params["GetInChI"] = "C" + cas_id
#             response = make_request_with_retry(nist_url, params)
#             if response is None:
#                 continue

#             num_created += 1
#             logging.info(
#                 f"Creating InChi key for id: {cas_id}. Total keys created {num_created}"
#             )
#             content = "{}\t{}\n".format(cas_id, response.content.decode("utf-8"))
#             file.write(content)


def scrap_inchi(cas_ls, params, data_dir):
    """Collect Inchi keys from NIST database and store them in txt format.

    Args:
        cas_ls: (list) CAS ids to download data for
        params: (dict) queries to be added to url
        data_dir: (string) path to store the data

    Returns:
        None
    """

    # Create file path for storing inchi keys
    inchi_path = os.path.join(data_dir, "inchi.txt")
    num_created = 0
    with open(inchi_path, "a") as file:
        content = "{}\t{}\n".format("cas_id", "inchi")
        file.write(content)

        for cas_id in cas_ls:
            params["GetInChI"] = "C" + cas_id
            response = requests.get(nist_url, params=params)

            num_created += 1
            logging.info(
                "Creating InChi key for id: {}. Total keys created {}".format(
                    cas_id, num_created
                )
            )
            content = "{}\t{}\n".format(cas_id, response.content.decode("utf-8"))
            file.write(content)


parser = argparse.ArgumentParser()
parser.add_argument(
    "--save_dir", default="./data", help="Directory path to store scrapped data"
)
parser.add_argument(
    "--cas_list",
    default="species.txt",
    help="File containing CAS number and formula of molecules",
)
parser.add_argument("--scrap_IR", default=True, help="Whether to download IR or not")
parser.add_argument("--scrap_MS", default=False, help="Whether to download MS or not")
parser.add_argument(
    "--scrap_InChi", default=True, help="Whether to download InChi or not"
)

args = parser.parse_args()

assert os.path.isfile(args.cas_list), "No file named {} exists".format(args.cas_list)

data_dir = args.save_dir
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

set_logger(data_dir, "scrap.log")

logging.info("Loading CAS file")
cas_df = pd.read_csv(
    args.cas_list, sep="\t", names=["name", "formula", "cas"], header=0
)
cas_df.dropna(subset=["cas"], inplace=True)
cas_df.cas = cas_df.cas.str.replace("-", "")

cas_ids = list(cas_df.cas)

logging.info("Scrap InChi keys")
if args.scrap_InChi:
    params = {}
    scrap_inchi(cas_ids, params, data_dir)

logging.info("Scrap IR spectra")
if args.scrap_IR:
    params = {"JCAMP": "", "Type": "IR", "Index": 0}
    scrap_data(cas_ids, params, data_dir)


logging.info("Scrap Mass spectra")
if args.scrap_MS:
    params = {"JCAMP": "", "Index": 0, "Type": "Mass"}
    scrap_data(cas_ids, params, data_dir)
