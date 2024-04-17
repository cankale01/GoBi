def process_file(file_path):
    sequence_ranges = {}

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            sequence_id = parts[0]
            start = int(parts[6])
            end = int(parts[7])

            if sequence_id not in sequence_ranges:
                sequence_ranges[sequence_id] = {'min': start, 'max': end}
            else:
                if start < sequence_ranges[sequence_id]['min']:
                    sequence_ranges[sequence_id]['min'] = start
                if end > sequence_ranges[sequence_id]['max']:
                    sequence_ranges[sequence_id]['max'] = end

    return sequence_ranges


def adjust_ranges(sequence_ranges):
    for sequence_id, values in sequence_ranges.items():
        values['min'] -= 1000
        values['max'] += 1000
    return sequence_ranges


def main():
    file_path = input("Enter the path to the input file: ")
    sequence_ranges = process_file(file_path)
    adjusted_ranges = adjust_ranges(sequence_ranges)

    for sequence_id, values in adjusted_ranges.items():
        print(f"Sequence ID: {sequence_id}, Smallest number: {values['min']}, Largest number: {values['max']}")


if __name__ == "__main__":
    main()
