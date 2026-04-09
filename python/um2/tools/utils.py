def string_to_lattice(text):
    rows = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([int(x) for x in line.split()])

    if not rows:
        raise ValueError("lattice string is empty")

    width = len(rows[0])
    for row in rows:
        if len(row) != width:
            raise ValueError("all lattice rows must have the same length")

    return rows