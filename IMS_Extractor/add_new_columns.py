from pathlib import Path
import yaml

new_extractor_columns: list


def load(data):
    data = Path(data)
    with data.open('r') as f:
        s = yaml.load(f, Loader=yaml.BaseLoader)
    global new_extractor_columns
    new_extractor_columns = s['new_extractor_columns']

    return s


def main():
    # Use to test loaded columns
    data_file = "new_extractor_columns.yaml"
    column_list = load(data_file)
    print(column_list)


if __name__ == "__main__":
    main()
