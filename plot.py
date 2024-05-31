import matplotlib.pyplot as plt
import sys

def parse_file(file_path):
    parsed_data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            parsed_data.append(tuple(map(float, parts)))
    return parsed_data


file_path = sys.argv[1]
data = parse_file(file_path)

shapes = []
xs = []
ys = []

for p in data:
    if p == ():
        shapes.append((xs, ys))
        xs = []
        ys = []
        continue

    x, y = p
    xs.append(x)
    ys.append(y)
shapes.append((xs, ys))

for xs, ys in shapes:
    if len(xs) == 1:
        plt.plot(xs, ys, marker="o")
    else:
        plt.plot(xs, ys)

plt.gca().set_aspect('equal')
plt.show()
