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

ss = []
xs = []
ys = []

for s, x, y in data:
    ss.append(s)
    xs.append(x)
    ys.append(y)

plt.plot(xs, ys)
plt.gca().set_aspect('equal')
plt.show()
