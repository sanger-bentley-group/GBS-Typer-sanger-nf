def get_items(input_file):

    list_of_items = []
    with open(input_file, "r") as file:
        for line in file:
            list_of_items.append(line.split('\n')[0])

    return list_of_items
