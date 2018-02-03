import csv

def parse_freesurfer_stats(filename, output):
    with open(filename) as f:
        lines = f.readlines()

    stats = [['StructureName', 'NumberVoxels', 'Volume', 'Mean', 'STD', 'Min', 'Max', 'Range']]
    index = 1
    for i in range(len(lines)):
        if lines[i][0] != "#" and lines[i][0] != "\n":
            line_stats = lines[i].split()
            stats.append([])
            stats[index].append(line_stats[4])
            stats[index].append(line_stats[2])
            stats[index].append(line_stats[3])
            stats[index].append(line_stats[5])
            stats[index].append(line_stats[6])
            stats[index].append(line_stats[7])
            stats[index].append(line_stats[8])
            stats[index].append(line_stats[9])
            index += 1


    with open(output, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(stats)
