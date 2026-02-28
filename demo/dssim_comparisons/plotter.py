#!/usr/bin/python3
import matplotlib.pyplot as plt
import sys

if len(sys.argv) <= 1:
    print("Please enter a path to a folder containing a DSSIM report.txt")
    sys.exit(-1)

dir = sys.argv[1]
data=dict()
with open(dir + '/report.txt', 'r') as file:
    for line in file:
        # Process each line here
        (val,file) = line.strip().split()
        val = float(val)
        file = file.split('/')[1]
        if not file.endswith(".png"):
            print("unknown item:",file)
        elif file.startswith('temp_'):
            data[file[5:-4]] = val
        elif file.startswith('temp'):
            q = int(file[4:7])
            key = file[8:-4]
            if key not in data:
                data[key]=[]
            data[key].append((q,val))

for key, value in data.items():
    if type(value) is list:
        xs = [p[0] for p in value]
        ys = [p[1] for p in value]
        plt.plot(xs,ys,label=key)
    else:
        plt.scatter(101,value,label=key)


# Add labels and a title
plt.yscale('log') # Set the y-axis to a logarithmic scale
plt.xlabel('Encoder Q-Factor')
plt.ylabel('DSSIM')
plt.title(dir)
plt.legend()

# Display the plot
plt.savefig(dir + '/report.png')
plt.show()
