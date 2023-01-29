import argparse
import json
import os


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=True)
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))
    else:
        settings = [{}]
    return settings


def main():
    settings = parse_arguments_to_settings()
    project_dir = settings['project_dir']
    with open('big_clean.sh', 'w') as f:
        for sample in os.listdir(project_dir):
            d = {
                'alignment_dir': os.path.join(project_dir, sample),
                'sample': sample,
            }
            new_line = """
alignment_dir={alignment_dir}
sample={sample}
""".format(**d)
            f.write(new_line)
            f.write(open("clean_intermediate_files_template.sh").read())


if __name__ == '__main__':
    main()
