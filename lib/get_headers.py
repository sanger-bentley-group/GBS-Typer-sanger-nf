import json

def read_header_json(header_file):
    with open(header_file, 'r') as file:
        header_json = file.read()

    header_dict = json.loads(header_json)

    return header_dict
