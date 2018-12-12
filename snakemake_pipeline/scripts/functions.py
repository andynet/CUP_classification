def get_file_paths(wildcards):

    file = None
    if wildcards.type == 'normal':
        file = config['normal_file']
    elif wildcards.type == 'tumor':
        file = config['tumor_file']
    else:
        print('Incorrect file type. Possible values = normal | tumor')
        print(wildcards)
        exit()

    with open(file) as f:
        lines = f.readlines()

    result = [None, None]
    for line in lines:
        if line.startswith(wildcards.sample):
            result = [path.strip() for path in line.split('\t')[1:3]]

    return result


def get_samples(path, number):

    result = []

    with open(path) as f:
        lines = f.readlines()

    for i in range(number):
        result.append(lines[i].split('\t')[0])

    return result


def get_ref_base(path):

    import os

    return os.path.basename(path).replace('.fa', '')
